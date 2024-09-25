#!/usr/bin/env python

import numpy as np

rdist = np.array(
    [
        ["++", 0.1],
        ["+-", 0.1],
        ["+x", 0.1],
        ["-+", 0.1],
        ["--", 0.1],
        ["-x", 0.1],
        ["x+", 0.1],
        ["x-", 0.1],
        ["xx", 0.2],
    ]
)

rdist = [
    rdist[:, 0],
    rdist[:, 1].astype(float),
]

dist = {
    "+": lambda: 1 + np.random.exponential(scale=0.03),
    "-": lambda: max(0, 1 - np.random.exponential(scale=0.21)),
    "x": lambda: 0,
}

c = np.random.choice(rdist[0], p=rdist[1], size=(1000, 4))


def apply_dist(x: str, axis: int) -> float:
    return dist[x[axis]]()


f = np.vectorize(apply_dist, otypes=[np.float64])

tables = [f(c, axis=axis) for axis in range(2)]

for i, table in enumerate(tables):
    np.save(f"fitness_table_{i}", table)

print(tables)
