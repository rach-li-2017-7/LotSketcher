from dataclasses import dataclass, replace
from typing import Self

import numpy as np
from shapely import Point


@dataclass(frozen=True)
class PhysicsPoint:
    pos: np.ndarray
    vel: np.ndarray
    key: str
    fixed: bool = False


@dataclass
class PointSystem:
    points: list[PhysicsPoint]

    def step(self, delta_t: float, forces: list[np.ndarray]) -> Self:
        """
        dead simple euler integration for now
        """

        new_points = []
        for point, force in zip(self.points, forces):
            if point.fixed:
                new_point = replace(point)
            else:
                new_vel = point.vel + force * delta_t
                new_pos = point.pos + new_vel * delta_t
                new_point = replace(point, vel=new_vel, pos=new_pos)
            new_points.append(new_point)

        return PointSystem(new_points)
