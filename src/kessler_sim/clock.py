"""Simulation clock: sim time in seconds since the scenario epoch."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Clock:
    t: float = 0.0              # sim seconds since epoch
    speed: float = 60.0         # sim seconds per wall-clock second (negative = backwards)
    running: bool = False
    min_t: float = 0.0          # the scenario epoch is the earliest time

    def play(self):
        self.running = True

    def pause(self):
        self.running = False

    def set_speed(self, speed: float):
        self.speed = float(speed)

    def seek(self, t: float) -> float:
        self.t = max(self.min_t, float(t))
        return self.t

    def advance(self, wall_dt: float) -> float:
        """Advance by ``wall_dt`` real seconds (no-op when paused). Returns the new sim time."""
        if self.running and wall_dt > 0:
            self.seek(self.t + self.speed * wall_dt)
            if self.t == self.min_t and self.speed < 0:
                self.running = False
        return self.t

    def status(self) -> dict:
        return {"t": self.t, "speed": self.speed, "running": self.running}
