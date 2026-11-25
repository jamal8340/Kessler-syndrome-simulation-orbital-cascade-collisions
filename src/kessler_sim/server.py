"""WebSocket server: connects the simulation engine to Unity (see TRANSPORT_ARCHITECTURE.md).

    uv run python -m kessler_sim.server --scenario catalogue

* A streaming loop advances the clock and sends ~20 binary state messages per second
  to every client (only when the time or the catalogue changed).
* Commands are JSON; every reply carries the command's ``request_id``.
* Slow requests (close-approach search, scenario load) run in a background thread, so
  state messages keep flowing. Commands that change the engine wait until no
  background job is reading it.
* Per client, JSON messages are always sent before the next state message, and a slow
  client only ever gets the newest state (older ones are dropped).
"""

from __future__ import annotations

import argparse
import asyncio
import json
import logging
import time
from collections import deque

from websockets.asyncio.server import serve
from websockets.exceptions import ConnectionClosed

from . import protocol
from .breakup import BreakupParams
from .engine import Simulation
from .scenario import available_scenarios, load_scenario

log = logging.getLogger("kessler_sim.server")

DEFAULT_PORT = 8765
DEFAULT_RATE_HZ = 20.0
MAX_SEARCH_WINDOW_S = 6 * 3600.0


class CommandError(Exception):
    """A command that cannot be executed; reported to the client as an ``error``."""


class Client:
    """One connection: an ordered JSON queue plus a 'latest state' slot."""

    def __init__(self, ws):
        self.ws = ws
        self.queue: deque[str] = deque()
        self.frame: bytes | None = None
        self.wake = asyncio.Event()

    def send_json(self, msg: dict):
        self.queue.append(protocol.dumps(msg))
        self.wake.set()

    def send_frame(self, data: bytes):
        self.frame = data              # replaces an unsent older frame
        self.wake.set()

    async def run(self):
        while True:
            await self.wake.wait()
            self.wake.clear()
            while self.queue:
                await self.ws.send(self.queue.popleft())
            if self.frame is not None:
                data, self.frame = self.frame, None
                await self.ws.send(data)


class SimulationServer:
    def __init__(self, sim: Simulation, host: str = "127.0.0.1", port: int = DEFAULT_PORT,
                 rate_hz: float = DEFAULT_RATE_HZ):
        self.sim = sim
        self.host, self.port = host, port
        self.period = 1.0 / rate_hz
        self.clients: set[Client] = set()
        self.jobs = 0
        self.idle: asyncio.Condition | None = None
        self._last_frame_key = None
        self._server = None
        self._stream_task = None

    # -- lifecycle ---------------------------------------------------------------

    async def start(self):
        self.idle = asyncio.Condition()
        self._server = await serve(self._handle, self.host, self.port, max_size=2 ** 22)
        self.port = self._server.sockets[0].getsockname()[1]
        self._stream_task = asyncio.create_task(self._stream())
        log.info("listening on ws://%s:%d (scenario '%s', %d objects)",
                 self.host, self.port, self.sim.name, len(self.sim.table))
        return self

    async def close(self):
        if self._stream_task:
            self._stream_task.cancel()
        if self._server:
            self._server.close()
            await self._server.wait_closed()

    async def serve_forever(self):
        await self.start()
        try:
            await asyncio.Future()
        finally:
            await self.close()

    # -- broadcasting ------------------------------------------------------------

    def broadcast(self, msg: dict, requester: Client | None = None):
        """Send to every client; only the requester sees the message's request_id."""
        anonymous = {k: v for k, v in msg.items() if k != "request_id"}
        for c in list(self.clients):
            c.send_json(msg if c is requester else anonymous)

    def broadcast_state(self, force: bool = False):
        key = (id(self.sim), self.sim.clock.t, self.sim.table.version)
        if not self.clients or (key == self._last_frame_key and not force):
            return
        data = protocol.state_message(self.sim)
        for c in list(self.clients):
            c.send_frame(data)
        self._last_frame_key = key

    def broadcast_changes(self, event=None, request_id=None, requester: Client | None = None):
        """events_changed -> breakup_event -> catalogue_update, in that order.

        Only the breakup_event is the reply (carries the request_id, for the requester).
        """
        changes = self.sim.take_changes()
        if changes.removed_events:
            self.broadcast(protocol.events_changed(changes.removed_events))
        if event is not None:
            self.broadcast(protocol.breakup_event(event, request_id), requester)
        if changes.removed or changes.updated or changes.added:
            self.broadcast(protocol.catalogue_update(self.sim, changes))
        self.broadcast_state(force=True)

    async def _stream(self):
        last = time.monotonic()
        next_tick = last
        while True:
            next_tick += self.period                       # fixed cadence, work time included
            await asyncio.sleep(max(0.0, next_tick - time.monotonic()))
            now = time.monotonic()
            if now - next_tick > self.period:              # fell behind: do not try to catch up
                next_tick = now
            try:
                clock = self.sim.clock
                was_running = clock.running
                clock.advance(now - last)
                if was_running and not clock.running:      # stopped at the epoch
                    self.broadcast(protocol.clock_status(self.sim))
                self.broadcast_state()
            except Exception:                               # never let the stream die
                log.exception("streaming step failed")
            last = now

    # -- engine access -----------------------------------------------------------

    async def _background(self, fn, *args):
        """Run a read-only job in a thread; engine changes wait until it is done."""
        async with self.idle:
            self.jobs += 1
        try:
            return await asyncio.to_thread(fn, *args)
        finally:
            async with self.idle:
                self.jobs -= 1
                self.idle.notify_all()

    # -- connections -------------------------------------------------------------

    async def _handle(self, ws):
        client = Client(ws)
        self.clients.add(client)
        sender = asyncio.create_task(client.run())
        log.info("client connected (%d total)", len(self.clients))
        try:
            async for raw in ws:
                if isinstance(raw, bytes):
                    client.send_json(protocol.error("binary messages are not accepted"))
                    continue
                try:
                    msg = json.loads(raw)
                    if not isinstance(msg, dict):
                        raise ValueError
                except ValueError:
                    client.send_json(protocol.error("message is not a JSON object"))
                    continue
                if msg.get("type") in ("list_close_approaches", "load_scenario"):
                    asyncio.create_task(self._dispatch(client, msg))     # slow: do not block
                else:
                    await self._dispatch(client, msg)
        except ConnectionClosed:
            pass
        finally:
            self.clients.discard(client)
            sender.cancel()
            log.info("client disconnected (%d left)", len(self.clients))

    async def _dispatch(self, client: Client, msg: dict):
        kind = msg.get("type")
        request_id = msg.get("request_id")
        handler = getattr(self, f"_cmd_{kind}", None) if isinstance(kind, str) else None
        if handler is None:
            client.send_json(protocol.error(f"unknown message type '{kind}'", request_id))
            return
        try:
            await handler(client, msg, request_id)
        except (CommandError, KeyError, ValueError, TypeError) as exc:
            reason = exc.args[0] if exc.args else type(exc).__name__
            client.send_json(protocol.error(str(reason), request_id))
        except ConnectionClosed:
            pass
        except Exception as exc:
            log.exception("command %s failed", kind)
            client.send_json(protocol.error(f"internal error: {exc}", request_id))

    # -- commands ----------------------------------------------------------------

    async def _cmd_hello(self, client, msg, request_id):
        client.send_json(protocol.session(self.sim, "hello", request_id))
        client.send_frame(protocol.state_message(self.sim))

    async def _cmd_clock(self, client, msg, request_id):
        clock = self.sim.clock
        action = msg.get("action")
        if action == "play":
            clock.play()
        elif action == "pause":
            clock.pause()
        elif action == "toggle":
            clock.pause() if clock.running else clock.play()
        elif action == "set_speed":
            clock.set_speed(float(msg["speed"]))
        else:
            raise CommandError("clock action must be play, pause, toggle or set_speed")
        self.broadcast(protocol.clock_status(self.sim, request_id), client)

    async def _cmd_seek(self, client, msg, request_id):
        self.sim.clock.seek(float(msg["t"]))
        self.broadcast(protocol.clock_status(self.sim, request_id), client)
        self.broadcast_state(force=True)

    async def _cmd_object_info(self, client, msg, request_id):
        t = msg.get("t")
        info = self.sim.object_info(int(msg["id"]), None if t is None else float(t))
        client.send_json(protocol.message("object_info", request_id, **info))

    async def _cmd_list_close_approaches(self, client, msg, request_id):
        sim = self.sim
        t_start = float(msg["t_start"]) if msg.get("t_start") is not None else sim.clock.t
        duration = float(msg.get("duration_s", 600.0))
        if not 0 < duration <= MAX_SEARCH_WINDOW_S:
            raise CommandError(f"duration_s must be in (0, {MAX_SEARCH_WINDOW_S:.0f}]")
        top = int(msg.get("top", 50))
        found = await self._background(sim.find_close_approaches, t_start, duration)
        client.send_json(protocol.message("close_approaches", request_id, t_start=t_start, duration_s=duration,
                                          total=len(found), approaches=found[:top]))

    async def _cmd_trigger_breakup(self, client, msg, request_id):
        params = BreakupParams(**(msg.get("params") or {}))
        t = msg.get("t")
        async with self.idle:
            await self.idle.wait_for(lambda: self.jobs == 0)
            try:
                event = self.sim.trigger_breakup(int(msg["id_a"]), int(msg["id_b"]),
                                                 None if t is None else float(t), params)
            except Exception:
                self.broadcast_changes()             # history may already have been rewritten
                raise
            self.broadcast_changes(event, request_id, client)

    async def _cmd_load_scenario(self, client, msg, request_id):
        name = str(msg["name"])
        if name not in available_scenarios():
            raise CommandError(f"unknown scenario '{name}'; available: {available_scenarios()}")
        new_sim = await asyncio.to_thread(load_scenario, name)
        new_sim.take_changes()                     # the session message carries everything
        async with self.idle:
            await self.idle.wait_for(lambda: self.jobs == 0)
            self.sim = new_sim
            self.broadcast(protocol.session(new_sim, "scenario_loaded", request_id), client)
            self.broadcast_state(force=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Kessler simulation WebSocket server")
    ap.add_argument("--scenario", default="catalogue", help=f"one of {available_scenarios()}")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=DEFAULT_PORT)
    ap.add_argument("--rate", type=float, default=DEFAULT_RATE_HZ, help="state messages per second")
    ap.add_argument("--speed", type=float, default=60.0, help="initial sim seconds per real second")
    ap.add_argument("--play", action="store_true", help="start with the clock running")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    log.info("loading scenario '%s'...", args.scenario)
    sim = load_scenario(args.scenario)
    sim.take_changes()
    sim.clock.set_speed(args.speed)
    if args.play:
        sim.clock.play()
    try:
        asyncio.run(SimulationServer(sim, args.host, args.port, args.rate).serve_forever())
    except KeyboardInterrupt:
        log.info("stopped")


if __name__ == "__main__":
    main()
