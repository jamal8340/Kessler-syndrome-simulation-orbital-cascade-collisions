"""End-to-end tests: a real server on a free port, talked to over WebSocket like Unity does."""

import asyncio
import json

import numpy as np
import pytest

from kessler_sim import Simulation, load_scenario
from kessler_sim.client import BackendConnection
from kessler_sim.server import SimulationServer
from test_protocol import three_body_sim
from test_sim_engine import T_HIT


def run_with_server(sim, scenario_fn, rate_hz=50.0, clients=1):
    async def main():
        server = await SimulationServer(sim, port=0, rate_hz=rate_hz).start()
        url = f"ws://127.0.0.1:{server.port}"
        try:
            conns = [BackendConnection(url) for _ in range(clients)]
            for c in conns:
                await c.__aenter__()
            try:
                return await scenario_fn(server, *conns)
            finally:
                for c in conns:
                    await c.__aexit__(None, None, None)
        finally:
            await server.close()
    return asyncio.run(main())


def alive(state):
    return int(np.sum(~np.isnan(state["positions"][:, 0])))


def test_hello_returns_session_and_first_state():
    async def scenario(server, conn):
        hello = await conn.request("hello")
        state = await conn.wait_for_state()
        return hello, state
    hello, state = run_with_server(load_scenario("iridium_cosmos_2009"), scenario)
    assert hello["type"] == "hello" and hello["request_id"] == 1
    assert hello["scenario"]["name"] == "iridium_cosmos_2009"
    assert len(hello["events"]) == 1
    assert len(state["positions"]) == len(hello["catalogue"])
    assert state["catalogue_version"] == hello["catalogue_version"]
    assert state["t"] == 0.0 and alive(state) == 2           # before the collision: the two satellites


def test_clock_play_streams_and_pause_stops():
    async def scenario(server, conn):
        await conn.request("hello")
        await conn.request("clock", action="set_speed", speed=100.0)
        status = await conn.request("clock", action="play")
        await asyncio.sleep(0.5)
        await conn.request("clock", action="pause")
        await asyncio.sleep(0.1)
        n_paused = len(conn.states)
        await asyncio.sleep(0.3)
        return status, [s["t"] for s in conn.states], n_paused, len(conn.states)
    status, times, n_paused, n_after = run_with_server(three_body_sim(), scenario, rate_hz=50.0)
    assert status["type"] == "clock_status" and status["running"] and status["speed"] == 100.0
    assert len(times) >= 15                                   # ~25 expected in 0.5 s at 50 Hz
    assert times == sorted(times) and times[-1] == pytest.approx(50.0, abs=15.0)
    assert n_after == n_paused                                # nothing new is sent while paused


def test_seek_sends_state_at_that_time():
    async def scenario(server, conn):
        await conn.request("hello")
        status = await conn.request("seek", t=1234.5)
        state = await conn.wait_for_state(lambda s: s["t"] == 1234.5)
        return status, state
    status, state = run_with_server(three_body_sim(), scenario)
    assert status["t"] == 1234.5 and not status["running"]
    assert alive(state) == 4


def test_breakup_messages_arrive_in_order_and_state_matches():
    async def scenario(server, conn):
        await conn.request("hello")
        reply = await conn.request("trigger_breakup", id_a=1, id_b=2, t=T_HIT)
        update = await conn.wait_for_message("catalogue_update")
        state = await conn.wait_for_state(lambda s: s["catalogue_version"] == update["catalogue_version"])
        return reply, update, state, [m["type"] for m in conn.messages]
    reply, update, state, order = run_with_server(three_body_sim(), scenario)
    assert reply["type"] == "breakup_event" and reply["request_id"] == 2
    assert order.index("breakup_event") < order.index("catalogue_update")
    assert len(update["added"]) == reply["n_fragments"]
    assert [e["id"] for e in update["updated"]] == [1, 2]
    assert all(e["t_destroyed"] == T_HIT for e in update["updated"])
    assert len(state["positions"]) == update["count"] == 4 + reply["n_fragments"]


def test_rewriting_history_over_the_network_keeps_client_in_sync():
    async def scenario(server, conn):
        hello = await conn.request("hello")
        await conn.request("trigger_breakup", id_a=1, id_b=2, t=T_HIT)
        await conn.request("trigger_breakup", id_a=1, id_b=3, t=T_HIT - 50)    # erases the first one
        await conn.wait_for_message("catalogue_update", count=2)
        again = await conn.request("hello")
        return hello, conn.messages, again
    hello, messages, again = run_with_server(three_body_sim(), scenario)

    types = [m["type"] for m in messages]
    second_reply = [i for i, m in enumerate(messages) if m["type"] == "breakup_event"][1]
    changed = types.index("events_changed")
    assert changed < second_reply                              # events_changed before the new event
    assert messages[changed]["removed"] == [0] and "request_id" not in messages[changed]

    # apply every update the way Unity must, starting from the first hello
    entries = hello["catalogue"]
    for m in messages:
        if m["type"] == "catalogue_update":
            removed = set(m["removed"])
            entries = [e for e in entries if e["id"] not in removed]
            upd = {e["id"]: e for e in m["updated"]}
            entries = [upd.get(e["id"], e) for e in entries] + m["added"]
    assert entries == again["catalogue"]


def test_close_approach_search_does_not_block_streaming():
    async def scenario(server, conn):
        await conn.request("hello")
        await conn.request("clock", action="play")
        before = len(conn.states)
        found = await conn.request("list_close_approaches", duration_s=600, top=5)
        return found, len(conn.states) - before
    found, frames_during_search = run_with_server(Simulation.from_catalogue(), scenario, rate_hz=20.0)
    assert found["type"] == "close_approaches" and found["total"] > 100
    assert len(found["approaches"]) == 5
    assert frames_during_search >= 10                         # search takes ~2 s; frames kept flowing


def test_errors_are_reported_and_server_keeps_working():
    async def scenario(server, conn):
        await conn.ws.send("not json")
        await conn.ws.send(json.dumps([1, 2]))
        await conn.ws.send(b"\x00\x01")
        replies = [await conn.request(kind, **fields) for kind, fields in [
            ("no_such_command", {}),
            ("trigger_breakup", {"id_a": 1, "id_b": 1}),
            ("trigger_breakup", {"id_a": 1, "id_b": 999}),
            ("trigger_breakup", {"id_a": 1, "id_b": 2, "params": {"bogus": 1}}),
            ("clock", {"action": "fly"}),
            ("list_close_approaches", {"duration_s": -5}),
            ("load_scenario", {"name": "nope"}),
            ("object_info", {"id": 12345}),
        ]]
        ok = await conn.request("hello")
        return replies, ok, [m for m in conn.messages if "request_id" not in m]
    replies, ok, unsolicited = run_with_server(three_body_sim(), scenario)
    assert all(r["type"] == "error" and r["reason"] for r in replies)
    assert [m["type"] for m in unsolicited[:3]] == ["error"] * 3
    assert ok["type"] == "hello"


def test_object_info_has_null_position_when_object_does_not_exist():
    async def scenario(server, conn):
        await conn.request("trigger_breakup", id_a=1, id_b=2, t=T_HIT)
        before = await conn.request("object_info", id=1, t=T_HIT - 1)
        after = await conn.request("object_info", id=1, t=T_HIT + 1)
        return before, after
    before, after = run_with_server(three_body_sim(), scenario)
    assert before["alive"] and before["position_km"] is not None
    assert not after["alive"] and after["position_km"] is None and after["t_destroyed"] == T_HIT


def test_two_clients_see_the_same_simulation_but_only_requester_gets_request_id():
    async def scenario(server, a, b):
        await a.request("hello")
        await b.request("hello")
        await a.request("trigger_breakup", id_a=1, id_b=2, t=T_HIT)
        event_b = await b.wait_for_message("breakup_event")
        await b.wait_for_message("catalogue_update")
        return event_b
    event_b = run_with_server(three_body_sim(), scenario, clients=2)
    assert "request_id" not in event_b and event_b["n_fragments"] > 0


def test_load_scenario_resets_everyone():
    async def scenario(server, conn):
        await conn.request("hello")
        loaded = await conn.request("load_scenario", name="iridium_cosmos_2009")
        state = await conn.wait_for_state(lambda s: len(s["positions"]) == len(loaded["catalogue"]))
        return loaded, state
    loaded, state = run_with_server(three_body_sim(), scenario)
    assert loaded["type"] == "scenario_loaded" and loaded["scenario"]["name"] == "iridium_cosmos_2009"
    assert len(loaded["events"]) == 1 and state["catalogue_version"] == loaded["catalogue_version"]
