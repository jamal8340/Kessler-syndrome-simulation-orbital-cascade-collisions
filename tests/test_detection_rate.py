import math, random, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src" / "collision_detection"))
import close_approach as ca

def make_pair(rng, window_s, d_miss, v_rel):
    # losowy kierunek ruchu względnego i losowe miejsce w LEO (różne pozycje względem siatki)
    u = [rng.gauss(0, 1) for _ in range(3)]
    n = math.sqrt(sum(c * c for c in u)); u = [c / n for c in u]
    w = [rng.gauss(0, 1) for _ in range(3)]
    dotuw = sum(a * b for a, b in zip(u, w))
    w = [b - dotuw * a for a, b in zip(u, w)]          # wektor prostopadły do u
    nw = math.sqrt(sum(c * c for c in w)); w = [c / nw for c in w]
    p = [rng.uniform(-7000, 7000) for _ in range(3)]
    t_c = rng.uniform(0, window_s)
    vA = [0.5 * v_rel * c for c in u]
    vB = [-0.5 * v_rel * c for c in u]
    rA_tca = p
    rB_tca = [p[i] + d_miss * w[i] for i in range(3)]
    rA0 = [rA_tca[i] - vA[i] * t_c for i in range(3)]
    rB0 = [rB_tca[i] - vB[i] * t_c for i in range(3)]
    return (tuple(rA0), tuple(vA)), (tuple(rB0), tuple(vB))

def detected(A, B, step_s, cell_km, horizon_s, threshold_km):
    (rA, vA), (rB, vB) = A, B
    t = 0.0
    while t < horizon_s:
        a = ca.add_scaled(rA, vA, t); b = ca.add_scaled(rB, vB, t)
        objs = [{"r": a}, {"r": b}]
        if ca.screen_candidate_pairs(objs, cell_km):      # czy siatka w ogóle zgłosi parę
            _, d_min, rel, _ = ca.compute_tca_and_miss_distance(a, vA, b, vB, step_s)
            if d_min <= threshold_km and rel >= ca.MIN_RELATIVE_SPEED_KMS:
                return True
        t += step_s
    return False

def run(step_s, cell_km, n=400, horizon_s=600.0, seed=1):
    rng = random.Random(seed)
    bins = [(0.5, 3), (3, 8), (8, 15)]
    for lo, hi in bins:
        hits = 0
        for _ in range(n):
            v = rng.uniform(lo, hi); d = rng.uniform(0, 5)
            A, B = make_pair(rng, horizon_s - step_s, d, v)
            hits += detected(A, B, step_s, cell_km, horizon_s, ca.DEFAULT_SCREENING_THRESHOLD_KM)
        print(f"step={step_s:>4}s cell={cell_km:>5}km v_rel {lo}-{hi} km/s: {100*hits/n:5.1f}%")

if __name__ == "__main__":
    print("Obecne ustawienia:"); run(ca.PROPAGATION_STEP_SECONDS, ca.CELL_SIZE_KM)
    print("Po poprawce:");       run(5.0, 100.0)
