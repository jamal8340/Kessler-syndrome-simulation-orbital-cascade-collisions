"""Command-line test client for the server (no Unity needed).

    uv run python -m kessler_sim.server --scenario catalogue      # in one terminal
    uv run python -m kessler_sim.client info                      # in another
    uv run python -m kessler_sim.client watch --seconds 5
    uv run python -m kessler_sim.client demo

``BackendConnection`` is also used by the tests: it talks to the server exactly
like Unity does.
"""

from __future__ import annotations

import argparse
import asyncio
import itertools
import json

import numpy as np
from websockets.asyncio.client import connect

from . import protocol
from .server import DEFAULT_PORT


class BackendConnection:
    def __init__(self, url: str):
        self.url = url
        self.ws = None
        self.messages: list[dict] = []          # every JSON message received, in order
        self.states: list[dict] = []            # every decoded state message, in order
        self._ids = itertools.count(1)
        self._waiters: dict[int, asyncio.Future] = {}
        self._new_message = asyncio.Event()
        self._reader = None

    async def __aenter__(self):
        self.ws = await connect(self.url, max_size=2 ** 26)
        self._reader = asyncio.create_task(self._read())
        return self

    async def __aexit__(self, *exc):
        self._reader.cancel()
        await self.ws.close()

    async def _read(self):
        async for raw in self.ws:
            if isinstance(raw, bytes):
                self.states.append(protocol.unpack_state(raw))
            else:
                msg = json.loads(raw)
                self.messages.append(msg)
                fut = self._waiters.pop(msg.get("request_id"), None)
                if fut is not None and not fut.done():
                    fut.set_result(msg)
            self._new_message.set()

    async def request(self, kind: str, timeout: float = 30.0, **fields) -> dict:
        """Send a command and wait for the reply carrying its request_id."""
        request_id = next(self._ids)
        fut = asyncio.get_running_loop().create_future()
        self._waiters[request_id] = fut
        await self.ws.send(json.dumps({"type": kind, "request_id": request_id, **fields}))
        return await asyncio.wait_for(fut, timeout)

    async def wait_for_message(self, kind: str, count: int = 1, timeout: float = 10.0) -> dict:
        """The ``count``-th JSON message of the given type (waits until it has arrived)."""
        async def _wait():
            while True:
                found = [m for m in self.messages if m["type"] == kind]
                if len(found) >= count:
                    return found[count - 1]
                self._new_message.clear()
                await self._new_message.wait()
        return await asyncio.wait_for(_wait(), timeout)

    async def wait_for_state(self, predicate=lambda s: True, timeout: float = 10.0) -> dict:
        async def _wait():
            seen = 0
            while True:
                for s in self.states[seen:]:
                    if predicate(s):
                        return s
                seen = len(self.states)
                self._new_message.clear()
                await self._new_message.wait()
        return await asyncio.wait_for(_wait(), timeout)


def _alive(state) -> int:
    return int(np.sum(~np.isnan(state["positions"][:, 0])))


async def cmd_info(url):
    async with BackendConnection(url) as conn:
        hello = await conn.request("hello")
        cats = {}
        for e in hello["catalogue"]:
            cats[e["category"]] = cats.get(e["category"], 0) + 1
        print(f"scenario: {hello['scenario']['name']} - {hello['scenario']['description']}")
        print(f"epoch: {hello['epoch_utc']}  frame: {hello['frame']}  catalogue version: {hello['catalogue_version']}")
        print(f"clock: {hello['clock']}")
        print(f"catalogue: {len(hello['catalogue'])} objects {cats}")
        print(f"events: {len(hello['events'])}")
        state = await conn.wait_for_state()
        print(f"first state: t = {state['t']:.1f} s, {_alive(state)} objects exist, "
              f"Earth angle {np.degrees(state['earth_rotation_rad']):.2f} deg, {len(state['positions'])} slots")


async def cmd_watch(url, seconds):
    async with BackendConnection(url) as conn:
        await conn.request("hello")
        await conn.request("clock", action="play")
        await asyncio.sleep(seconds)
        await conn.request("clock", action="pause")
        ts = [s["t"] for s in conn.states]
        print(f"received {len(conn.states)} state messages in {seconds:g} s "
              f"({len(conn.states) / seconds:.1f} per second); sim time {ts[0]:.1f} -> {ts[-1]:.1f} s")


async def cmd_demo(url):
    async with BackendConnection(url) as conn:
        hello = await conn.request("hello")
        print(f"connected: scenario '{hello['scenario']['name']}', {len(hello['catalogue'])} objects")
        print("searching close approaches in the next 10 minutes...")
        found = await conn.request("list_close_approaches", duration_s=600, top=3)
        if not found["approaches"]:
            print("none found")
            return
        top = found["approaches"][0]
        print(f"closest: {top['name_1']} x {top['name_2']}, {top['miss_distance_km']:.3f} km at t = {top['tca_s']:.1f} s")
        reply = await conn.request("trigger_breakup", id_a=top["id_1"], id_b=top["id_2"], t=top["tca_s"])
        print(f"breakup #{reply['event_id']}: {reply['n_fragments']} fragments "
              f"({reply['fragments_per_parent']}), catastrophic = {reply['catastrophic']}")
        update = await conn.wait_for_message("catalogue_update")
        print(f"catalogue_update: +{len(update['added'])} added, {len(update['updated'])} updated, "
              f"{len(update['removed'])} removed -> version {update['catalogue_version']}")
        for t in (top["tca_s"] - 60, top["tca_s"] + 60):
            await conn.request("seek", t=t)
            s = await conn.wait_for_state(lambda s, t=t: abs(s["t"] - t) < 1e-6
                                          and s["catalogue_version"] == update["catalogue_version"])
            print(f"t = {t:8.1f} s: {_alive(s)} objects exist (of {len(s['positions'])} slots)")
        info = await conn.request("object_info", id=top["id_1"], t=top["tca_s"] + 60)
        print(f"object_info of {info['name']} after the breakup: alive = {info['alive']}, "
              f"position = {info['position_km']}, destroyed at t = {info['t_destroyed']:.1f} s")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Test client for the Kessler simulation server")
    ap.add_argument("--url", default=f"ws://127.0.0.1:{DEFAULT_PORT}")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("info", help="print the session info and the first state message")
    watch = sub.add_parser("watch", help="play the clock and count state messages")
    watch.add_argument("--seconds", type=float, default=5.0)
    sub.add_parser("demo", help="find the closest approach, collide it, seek across it")
    args = ap.parse_args(argv)
    if args.cmd == "info":
        asyncio.run(cmd_info(args.url))
    elif args.cmd == "watch":
        asyncio.run(cmd_watch(args.url, args.seconds))
    else:
        asyncio.run(cmd_demo(args.url))


if __name__ == "__main__":
    main()
