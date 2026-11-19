"""
Penyelesaian: model SPL A.T = b dan tiga metode (Gauss-Seidel, Jacobi, Eliminasi Gauss).
Tidak ada tampilan di berkas ini; semua hasil dikembalikan sebagai data.

Model : suhu tiap titik dalam = rata-rata 4 tetangganya  ->  4T - (jumlah tetangga) = 0  ->  A.T = b
Nomor : T1..T(n*n) berurutan per baris dari kiri-atas (T1 = pojok kiri-atas, T2 di kanannya, T4 di bawahnya).
"""
import time
import numpy as np

from konfigurasi import TOP, SIDE, BOTTOM, TOL, MAXIT, ITER_SINGKAT


# =====================================================================================
# MODEL: suhu tepi dan SPL  A.T = b
# =====================================================================================
def boundary(n, top=TOP, left=SIDE, right=SIDE, bottom=BOTTOM):
    """Empat daftar suhu tepi (atas, bawah, kiri, kanan). top/bottom boleh angka atau fungsi x (0..1)."""
    xs = [(c + 1) / (n + 1) for c in range(n)]
    f = lambda v: [v(x) for x in xs] if callable(v) else [float(v)] * n
    g = lambda v: [float(v)] * n
    return f(top), f(bottom), g(left), g(right)


def build_system(n, bnd):
    """Susun matriks A (N x N) dan vektor b untuk grid n x n titik dalam, N = n*n."""
    top, bottom, left, right = bnd
    N = n * n
    A = np.zeros((N, N))
    b = np.zeros(N)
    for r in range(n):
        for c in range(n):
            k = r * n + c
            A[k, k] = 4.0
            if r == 0:
                b[k] += top[c]
            else:
                A[k, k - n] = -1.0
            if r == n - 1:
                b[k] += bottom[c]
            else:
                A[k, k + n] = -1.0
            if c == 0:
                b[k] += left[r]
            else:
                A[k, k - 1] = -1.0
            if c == n - 1:
                b[k] += right[r]
            else:
                A[k, k + 1] = -1.0
    return A, b


def tetangga(A):
    """Daftar indeks (0-based) tetangga tiap titik, dibaca dari matriks A."""
    N = A.shape[0]
    return [[j for j in range(N) if j != i and A[i, j] != 0] for i in range(N)]


def padded_grid(n, bnd, init):
    """Grid (n+2) x (n+2): tepi diisi suhu tepi, titik dalam diisi tebakan awal."""
    top, bottom, left, right = bnd
    g = [[float(init)] * (n + 2) for _ in range(n + 2)]
    for c in range(n):
        g[0][c + 1] = top[c]
        g[n + 1][c + 1] = bottom[c]
    for r in range(n):
        g[r + 1][0] = left[r]
        g[r + 1][n + 1] = right[r]
    return g


def interior(g, n):
    return np.array([g[r][1:n + 1] for r in range(1, n + 1)])


# =====================================================================================
# METODE ITERATIF
# =====================================================================================
def gauss_seidel(n, bnd, init, tol=TOL, maxit=MAXIT, keep=False):
    """Update di tempat: nilai atas dan kiri sudah yang terbaru pada iterasi yang sama.
    Mengembalikan (T, jumlah iterasi, riwayat selisih terbesar, riwayat semua T bila keep).
    tol = 0 berarti tidak pernah berhenti karena selisih (berhenti hanya di maxit)."""
    g = padded_grid(n, bnd, init)
    hist, snaps = [], []
    it = 0
    while it < maxit:
        it += 1
        dmax = 0.0
        for r in range(1, n + 1):
            up, row, dn = g[r - 1], g[r], g[r + 1]
            for c in range(1, n + 1):
                new = (up[c] + dn[c] + row[c - 1] + row[c + 1]) * 0.25
                d = abs(new - row[c])
                if d > dmax:
                    dmax = d
                row[c] = new
        hist.append(dmax)
        if keep:
            snaps.append(interior(g, n).flatten())
        if dmax < tol:
            break
    return interior(g, n), it, hist, snaps


def jacobi(n, bnd, init, tol=TOL, maxit=MAXIT, keep=False):
    """Semua titik dihitung dari nilai iterasi sebelumnya (dua grid)."""
    old = padded_grid(n, bnd, init)
    new = [row[:] for row in old]
    hist, snaps = [], []
    it = 0
    while it < maxit:
        it += 1
        dmax = 0.0
        for r in range(1, n + 1):
            up, row, dn, out = old[r - 1], old[r], old[r + 1], new[r]
            for c in range(1, n + 1):
                v = (up[c] + dn[c] + row[c - 1] + row[c + 1]) * 0.25
                d = abs(v - row[c])
                if d > dmax:
                    dmax = d
                out[c] = v
        old, new = new, old
        hist.append(dmax)
        if keep:
            snaps.append(interior(old, n).flatten())
        if dmax < tol:
            break
    return interior(old, n), it, hist, snaps


# =====================================================================================
# ELIMINASI GAUSS
# =====================================================================================
def gauss_elimination(A, b):
    """Eliminasi Gauss (pivoting parsial) + substitusi mundur. Hanya mengembalikan hasil akhir."""
    A = A.copy()
    b = b.copy()
    N = len(b)
    for k in range(N - 1):
        p = k + int(np.argmax(np.abs(A[k:, k])))
        if A[p, k] == 0.0:
            raise ValueError("matriks singular")
        if p != k:
            A[[k, p]] = A[[p, k]]
            b[[k, p]] = b[[p, k]]
        m = A[k + 1:, k] / A[k, k]
        A[k + 1:, k:] -= np.outer(m, A[k, k:])
        b[k + 1:] -= m * b[k]
    x = np.zeros(N)
    for i in range(N - 1, -1, -1):
        x[i] = (b[i] - A[i, i + 1:] @ x[i + 1:]) / A[i, i]
    return x


def eliminasi_langkah(A, b):
    """Eliminasi Gauss yang mencatat setiap langkah (untuk grid kecil).
    Mengembalikan dict: awal = [A|b], langkah = daftar langkah, U, c, x."""
    N = len(b)
    M = np.hstack([A, b[:, None]]).astype(float)
    awal = M.copy()
    langkah = []
    for k in range(N - 1):
        p = k + int(np.argmax(np.abs(M[k:, k])))
        if M[p, k] == 0.0:
            raise ValueError("matriks singular")
        tukar = None
        if p != k:
            M[[k, p]] = M[[p, k]]
            tukar = (k, p)
        pivot = M[k, k]
        ops = []
        for i in range(k + 1, N):
            if abs(M[i, k]) > 1e-12:
                aik = M[i, k]
                m = aik / pivot
                ops.append((i, aik, m))
                M[i, :] -= m * M[k, :]
                M[i, k] = 0.0
        langkah.append({"k": k, "pivot": pivot, "tukar": tukar, "ops": ops, "M": M.copy()})
    U = M[:, :N].copy()
    c = M[:, N].copy()
    x = np.zeros(N)
    for i in range(N - 1, -1, -1):
        x[i] = (c[i] - U[i, i + 1:] @ x[i + 1:]) / U[i, i]
    return {"awal": awal, "langkah": langkah, "U": U, "c": c, "x": x}


# =====================================================================================
# ALAT BANTU
# =====================================================================================
def timed(fn, repeat=1):
    """Waktu terbaik dari beberapa kali jalan (detik) + hasil terakhir."""
    best, res = float("inf"), None
    for _ in range(repeat):
        t0 = time.perf_counter()
        res = fn()
        best = min(best, time.perf_counter() - t0)
    return best, res


def residual_inf(A, T, b):
    return float(np.max(np.abs(A @ np.asarray(T).flatten() - b)))


def jalankan(n, init, repeat=1, keep=False, tol=TOL, maxit=MAXIT):
    """Jalankan tiga metode pada grid n x n. Mengembalikan dict hasil."""
    bnd = boundary(n)
    A, b = build_system(n, bnd)
    t_el, x_el = timed(lambda: gauss_elimination(A, b), repeat)
    hasil = {"n": n, "N": n * n, "init": float(init), "A": A, "b": b, "ref": x_el, "singkat": tol == 0.0,
             "eliminasi": {"T": x_el, "waktu": t_el, "residual": residual_inf(A, x_el, b)}}
    for nama, fn in (("gs", gauss_seidel), ("jacobi", jacobi)):
        t, (T, it, hist, snaps) = timed(lambda: fn(n, bnd, init, tol=tol, maxit=maxit, keep=keep), repeat)
        hasil[nama] = {"T": T.flatten(), "iterasi": it, "riwayat": hist, "snaps": snaps, "waktu": t,
                       "error": float(np.max(np.abs(T.flatten() - x_el))), "residual": residual_inf(A, T, b)}
    return hasil


_cache = {}


def hitung_grid(n, init=0.0, singkat=False, pesan=True):
    """Hitung (atau ambil dari cache) tiga metode pada grid n x n.
    singkat=True : Gauss-Seidel dan Jacobi hanya ITER_SINGKAT iterasi (tidak sampai konvergen).
    singkat=False: sampai selisih terbesar < TOL."""
    key = (n, float(init), singkat)
    if key not in _cache:
        if pesan and n > 3:
            print(f"  Menghitung grid {n}x{n} ({n * n} persamaan), mohon tunggu...")
        _cache[key] = jalankan(n, init,
                               repeat=500 if n <= 3 else 1,
                               keep=(n <= 3 or singkat),
                               tol=0.0 if singkat else TOL,
                               maxit=ITER_SINGKAT if singkat else MAXIT)
    return _cache[key]
