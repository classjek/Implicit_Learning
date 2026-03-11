import sys
import numpy as np
import scipy.sparse as sp
from collections import defaultdict
import cvxpy as cp


def parse_sdpa(path):
    with open(path) as f:
        raw = [l.strip() for l in f
               if l.strip() and not l.strip().startswith('*')]

    idx = 0
    mDim      = int(raw[idx]);           idx += 1
    nBlock    = int(raw[idx]);           idx += 1
    bsizes    = list(map(int, raw[idx].split())); idx += 1

    # b-vector (objective): all on one line OR split across lines
    b_vals = []
    while len(b_vals) < mDim:
        b_vals.extend(map(float, raw[idx].split()))
        idx += 1
    b = np.array(b_vals[:mDim])

    # matrix entries
    entries = []
    for line in raw[idx:]:
        parts = line.split()
        if len(parts) == 5:
            entries.append((int(parts[0]), int(parts[1]),
                            int(parts[2]), int(parts[3]),
                            float(parts[4])))
    return mDim, nBlock, bsizes, b, entries


def solve(dat_file, solver=cp.MOSEK, verbose=True):
    print(f"Parsing {dat_file} ...")
    mDim, nBlock, bsizes, b_obj, entries = parse_sdpa(dat_file)
    print(f"  mDim={mDim}, nBlock={nBlock}, block sizes={bsizes}")
    print(f"  Total entries: {len(entries)}")

    y = cp.Variable(mDim, name="y")
    constraints = []

    # group entries by block number
    by_block = defaultdict(list)
    for (v, blk, r, c, val) in entries:
        by_block[blk].append((v, r, c, val))

    print("Building SDP constraints ...")
    for blk_idx, bsize in enumerate(bsizes):
        blk = blk_idx + 1
        n   = abs(bsize)
        is_lp = (bsize < 0)
        ents = by_block.get(blk, [])
        if not ents:
            continue

        # Build A_sparse (n² × mDim) and b0 (n²) such that
        #   vec( F_0 + Σ F_j y_j ) = A_sparse @ y + b0
        b0      = np.zeros(n * n)
        rows, cols, vals = [], [], []

        for (v, r, c, val) in ents:
            ri, ci = r - 1, c - 1
            if v == 0:                       # constant F_0
                b0[ri * n + ci] += val
                if ri != ci:
                    b0[ci * n + ri] += val
            else:                            # coefficient F_v
                rows.append(ri * n + ci);  cols.append(v - 1);  vals.append(val)
                if ri != ci:
                    rows.append(ci * n + ri); cols.append(v - 1); vals.append(val)

        A = sp.csr_matrix((vals, (rows, cols)), shape=(n * n, mDim))

        vec_expr = A @ y + b0                              # n² affine expression
        mat_expr = cp.reshape(vec_expr, (n, n))
        mat_sym  = (mat_expr + mat_expr.T) / 2            # ensure symmetry

        if is_lp:
            # each diagonal entry ≥ 0
            for i in range(n):
                d     = A.getrow(i * n + i)
                const = b0[i * n + i]
                if d.nnz > 0:
                    constraints.append(float(const) + cp.reshape(d @ y, ()) >= 0)
                # if nnz==0 and const < 0, log a warning
        else:
            constraints.append(mat_sym >> 0)   # PSD

    # E[1] = 1: the constant moment is always 1 (anchors the scale of the problem)
    constraints.append(y[0] == 1.0)
    prob = cp.Problem(cp.Minimize(0), constraints)  # pure feasibility

    print(f"Solving with {solver} ...")
    prob.solve(solver=solver, verbose=verbose)

    print(f"\nStatus:    {prob.status}")
    print(f"Objective: {prob.value}")
    return prob


if __name__ == "__main__":
    dat = (sys.argv[1] if len(sys.argv) > 1
           else "/workspace/Implicit_Learning/data/sparsepop_output_test.dat-s")
    solve(dat)