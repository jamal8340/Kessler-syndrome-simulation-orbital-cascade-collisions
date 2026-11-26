"""Packed piecewise orbital histories, with a vectorized latest-segment fast path."""
from __future__ import annotations

import numpy as np

from . import drag


class _Track:
    def __init__(self):
        self.data = np.empty((4, drag.WIDTH))
        self.count = 0

    def append(self, segment):
        if self.count == len(self.data):
            new = np.empty((self.count * 2, drag.WIDTH))
            new[:self.count] = self.data
            self.data = new
        self.data[self.count] = segment
        self.count += 1

    def at(self, t):
        idx = np.searchsorted(self.data[:self.count, drag.START], t, side="right") - 1
        return None if idx < 0 else self.data[idx]

    def pop_at(self, t):
        self.count = int(np.searchsorted(self.data[:self.count, drag.START], t, side="right"))
        if not self.count:
            return None
        self.count -= 1
        result = self.data[self.count].copy()
        if self.count < len(self.data) // 4:
            self.data = self.data[:max(4, self.count * 2)].copy()
        return result


class OrbitHistory:
    def __init__(self):
        self.latest = np.empty((0, drag.WIDTH))
        self._older: dict[int, _Track] = {}

    def add(self, n):
        self.latest = np.concatenate((self.latest, np.full((n, drag.WIDTH), np.nan)))

    def remove(self, ids, keep):
        for oid in ids[~keep]:
            self._older.pop(int(oid), None)
        self.latest = self.latest[keep]

    def append(self, ids, rows, segments):
        previous = self.latest[rows]
        for oid, segment in zip(ids[rows][np.isfinite(previous[:, drag.START])],
                                previous[np.isfinite(previous[:, drag.START])]):
            key = int(oid)
            track = self._older.get(key)
            if track is None:
                track = self._older[key] = _Track()
            track.append(segment)
        self.latest[rows] = segments

    def at(self, ids, rows, t):
        segments = self.latest[rows].copy()
        for j in np.flatnonzero(segments[:, drag.START] > t):
            track = self._older.get(int(ids[rows[j]]))
            segment = track.at(t) if track is not None else None
            segments[j] = np.nan if segment is None else segment
        return segments

    def truncate_after(self, ids, t):
        for row in np.flatnonzero(self.latest[:, drag.START] > t):
            track = self._older.get(int(ids[row]))
            segment = track.pop_at(t) if track is not None else None
            self.latest[row] = np.nan if segment is None else segment
            if track is not None and not track.count:
                del self._older[int(ids[row])]

    @property
    def nbytes(self):
        """Allocated NumPy buffers (not Python object overhead)."""
        return self.latest.nbytes + sum(track.data.nbytes for track in self._older.values())

    @property
    def segment_count(self):
        return int(np.isfinite(self.latest[:, drag.START]).sum()) + sum(t.count for t in self._older.values())
