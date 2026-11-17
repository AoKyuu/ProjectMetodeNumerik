import os
import time
from decimal import Decimal
import numpy as np

TOP, SIDE, BOTTOM = 85.0, 40.0, 25.0   # suhu sisi atas, kiri=kanan, bawah (derajat Celsius)
TOL = 1e-6
MAXIT = 200000
FOLDER_GAMBAR = "gambar_program"



# 1. MODEL: suhu tepi dan SPL  A.T = b

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



# 2. TIGA METODE

def gauss_seidel(n, bnd, init, tol=TOL, maxit=MAXIT, keep=False):
    """Update di tempat: nilai atas dan kiri sudah yang terbaru pada iterasi yang sama.
    Mengembalikan (T, jumlah iterasi, riwayat selisih terbesar, riwayat semua T bila keep)."""
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


def gauss_elimination(A, b):
    """Eliminasi Gauss dengan pivoting parsial + substitusi mundur."""
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



# 3. ALAT BANTU (waktu, residual, format angka berkoma)

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


def fx(x, d=4):
    """Angka dengan d angka di belakang koma, memakai koma desimal."""
    return f"{x:.{d}f}".replace(".", ",")


def sci(x):
    """Angka kecil ditulis biasa (tanpa huruf e), 2 angka penting, mis. 0,00000059."""
    return format(Decimal(f"{x:.1e}"), "f").replace(".", ",")


def ms(detik):
    """Waktu dalam satuan yang enak dibaca."""
    return f"{detik * 1e6:,.0f} mikrodetik".replace(",", ".") if detik < 1e-3 else f"{detik:.3f} detik".replace(".", ",")


def garis(ch="-", n=78):
    print(ch * n)


def tanya_tebakan():
    s = input("Tebakan awal semua titik dalam (derajat C) [Enter = 0]: ").strip().replace(",", ".")
    if s == "":
        return 0.0
    try:
        return float(s)
    except ValueError:
        print("  Input bukan angka, dipakai 0.")
        return 0.0



# 4. MENJALANKAN SATU GRID (dipakai menu 1, 2, dan 3)

def jalankan(n, init, repeat, keep=False):
    """Jalankan tiga metode pada grid n x n. Mengembalikan dict hasil."""
    bnd = boundary(n)
    A, b = build_system(n, bnd)
    t_el, x_el = timed(lambda: gauss_elimination(A, b), repeat)
    hasil = {"n": n, "N": n * n, "init": init, "A": A, "b": b, "ref": x_el,
             "eliminasi": {"T": x_el, "waktu": t_el, "residual": residual_inf(A, x_el, b)}}
    for nama, fn in (("gs", gauss_seidel), ("jacobi", jacobi)):
        t, (T, it, hist, snaps) = timed(lambda: fn(n, bnd, init, keep=keep), repeat)
        hasil[nama] = {"T": T.flatten(), "iterasi": it, "riwayat": hist, "snaps": snaps, "waktu": t,
                       "error": float(np.max(np.abs(T.flatten() - x_el))), "residual": residual_inf(A, T, b)}
    return hasil


# cache supaya grafik tidak menghitung ulang
_cache = {}


def ambil(n, init=0.0):
    key = (n, init)
    if key not in _cache:
        if n > 3:
            print(f"  Menghitung grid {n}x{n} ({n * n} persamaan), mohon tunggu...")
        _cache[key] = jalankan(n, init, repeat=500 if n <= 3 else 1, keep=(n == 3))
    return _cache[key]


def ringkasan(h):
    """Tabel ringkas jumlah iterasi, waktu, error, residual."""
    print(f"\n{'Metode':<17}{'Iterasi':>8}{'Waktu':>19}{'Error maks (C)':>19}{'Residual maks':>19}")
    garis("-", 82)
    for nama, kunci in (("Gauss-Seidel", "gs"), ("Jacobi", "jacobi")):
        r = h[kunci]
        print(f"{nama:<17}{r['iterasi']:>8}{ms(r['waktu']):>19}{sci(r['error']):>19}{sci(r['residual']):>19}")
    e = h["eliminasi"]
    print(f"{'Eliminasi Gauss':<17}{'langsung':>8}{ms(e['waktu']):>19}{'(acuan)':>19}{sci(e['residual']):>19}")
    garis("-", 82)
    print("Error = selisih mutlak terbesar terhadap Eliminasi Gauss.  Residual = nilai mutlak terbesar A.T - b.")
    print("Syarat berhenti: selisih terbesar antar-iterasi < 0,000001 (bukan jaminan error < 0,000001).")



# MENU 1: GRID 3x3

def tampil_sistem_3x3(A, b):
    print("\nPersamaan tiap titik dalam (T = rata-rata 4 tetangga, suhu tepi sudah dimasukkan):")
    for k in range(9):
        tetangga = [f"T{j + 1}" for j in range(9) if A[k, j] == -1.0]
        konst = int(b[k]) if b[k] == int(b[k]) else b[k]
        suku = ([str(konst)] if konst != 0 else []) + tetangga
        print(f"  T{k + 1} = ({' + '.join(suku)}) / 4")
    print("\nMatriks A (9x9):")
    for k in range(9):
        print("  [" + " ".join(f"{int(v):>2}" for v in A[k]) + " ]")
    print("Vektor b:", [int(v) for v in b])


def tabel_iterasi(h, kunci):
    r = h[kunci]
    nama = "Gauss-Seidel" if kunci == "gs" else "Jacobi"
    print(f"\nTabel iterasi {nama} (tebakan awal {fx(h['init'], 1)} C), 7 angka di belakang koma")
    kepala = "Iter " + "".join(f"{'T' + str(i + 1):>13}" for i in range(9)) + f"{'Selisih maks':>14}"
    print(kepala)
    garis("-", len(kepala))
    print(f"{0:<5}" + "".join(f"{fx(h['init'], 7):>13}" for _ in range(9)))
    for i, (snap, d) in enumerate(zip(r["snaps"], r["riwayat"]), start=1):
        print(f"{i:<5}" + "".join(f"{fx(v, 7):>13}" for v in snap) + f"{fx(d, 7):>14}")
    print(f"Berhenti pada iterasi {r['iterasi']}.")


def menu_3x3():
    garis("=")
    print("GRID 3x3  (9 titik dalam, 9 persamaan, matriks A 9x9)")
    garis("=")
    print("Suhu sisi: atas 85 C, kiri dan kanan 40 C, bawah 25 C")
    init = tanya_tebakan()
    h = ambil(3, init)
    tampil_sistem_3x3(h["A"], h["b"])
    print("\nHasil suhu stabil (derajat C):")
    kepala = f"{'Titik':<7}{'Gauss-Seidel':>16}{'Jacobi':>16}{'Eliminasi Gauss':>18}"
    print(kepala)
    garis()
    for i in range(9):
        print(f"T{i + 1:<6}{fx(h['gs']['T'][i], 4):>16}{fx(h['jacobi']['T'][i], 4):>16}{fx(h['eliminasi']['T'][i], 4):>18}")
    garis()
    ringkasan(h)
    while True:
        s = input("\nTampilkan tabel semua iterasi?  1 = Gauss-Seidel, 2 = Jacobi, 0 = kembali: ").strip()
        if s == "1":
            tabel_iterasi(h, "gs")
        elif s == "2":
            tabel_iterasi(h, "jacobi")
        elif s == "0" or s == "":
            return
        else:
            print("  Pilihan tidak dikenal.")



# MENU 2: GRID 39x39

def menu_39x39():
    garis("=")
    print("GRID 39x39  (1.521 titik dalam, 1.521 persamaan, matriks A 1.521 x 1.521)")
    garis("=")
    print("Suhu sisi: atas 85 C, kiri dan kanan 40 C, bawah 25 C")
    print("Hasil tidak ditabelkan (T1 sampai T1521 terlalu banyak); yang ditampilkan ringkasannya.")
    init = tanya_tebakan()
    h = ambil(39, init)
    ringkasan(h)
    T = h["eliminasi"]["T"].reshape(39, 39)
    print(f"\nSuhu titik tengah (baris 20, kolom 20): {fx(T[19, 19], 4)} C")
    print(f"Suhu terendah {fx(T.min(), 4)} C, tertinggi {fx(T.max(), 4)} C (di dalam rentang 25-85 C: "
          f"{'ya' if 25 <= T.min() and T.max() <= 85 else 'tidak'})")
    print("Untuk melihat gambarnya, pilih menu 3.")
    input("\nTekan Enter untuk kembali ke menu...")



# MENU 3: GRAFIK / ILUSTRASI

def siapkan_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    ink, ink2, grid = "#0b0b0b", "#52514e", "#e3e2dd"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": grid, "axes.labelcolor": ink2, "xtick.color": ink2,
                         "ytick.color": ink2, "text.color": ink, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titlesize": 11, "axes.titleweight": "bold", "axes.titlelocation": "left"})
    cmap = LinearSegmentedColormap.from_list("panas", ["#fff1e6", "#fdbb84", "#e6550d", "#7f2704"])
    return plt, cmap, ink, ink2, grid


def simpan(fig, nama):
    os.makedirs(FOLDER_GAMBAR, exist_ok=True)
    path = os.path.join(FOLDER_GAMBAR, nama)
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    print(f"  tersimpan: {os.path.abspath(path)}")


def padded_hasil(T, n):
    pad = np.full((n + 2, n + 2), np.nan)
    pad[1:-1, 1:-1] = T.reshape(n, n)
    pad[0, 1:-1], pad[-1, 1:-1], pad[1:-1, 0], pad[1:-1, -1] = TOP, BOTTOM, SIDE, SIDE
    return pad


def gambar_pelat_3x3():
    plt, cmap, ink, ink2, grid = siapkan_matplotlib()
    h = ambil(3)
    pad = padded_hasil(h["eliminasi"]["T"], 3)
    fig, ax = plt.subplots(figsize=(5.2, 4.6), dpi=200)
    im = ax.imshow(pad, cmap=cmap, vmin=25, vmax=85)
    for r in range(5):
        for c in range(5):
            if np.isnan(pad[r, c]):
                continue
            dalam = 1 <= r <= 3 and 1 <= c <= 3
            k = (r - 1) * 3 + c
            v = fx(pad[r, c], 4) if dalam else fx(pad[r, c], 0)
            ax.text(c, r, f"T{k}\n{v}" if dalam else v, ha="center", va="center", fontsize=9 if dalam else 9.5,
                    color="white" if (dalam and pad[r, c] > 55) else ink, fontweight="bold" if dalam else "normal")
    for i in range(6):
        ax.axhline(i - 0.5, color="white", lw=1.5)
        ax.axvline(i - 0.5, color="white", lw=1.5)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Suhu stabil contoh 3×3 (sel tepi = suhu sisi yang diketahui)", pad=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04); cb.set_label("Suhu (°C)"); cb.outline.set_visible(False)
    simpan(fig, "gambar1_pelat_3x3.png"); plt.close(fig)


def grafik_konvergensi(n, nama_file):
    plt, cmap, ink, ink2, grid = siapkan_matplotlib()
    c1, c2 = "#2a78d6", "#eb6834"      # biru = Gauss-Seidel, oranye = Jacobi
    h = ambil(n)
    hg, hj = h["gs"]["riwayat"], h["jacobi"]["riwayat"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6), dpi=200)
    kecil = n <= 3
    ax.plot(range(1, len(hj) + 1), hj, color=c2, lw=2, label="Jacobi", **({"marker": "o", "ms": 3.5} if kecil else {}))
    ax.plot(range(1, len(hg) + 1), hg, color=c1, lw=2, label="Gauss-Seidel", **({"marker": "o", "ms": 3.5} if kecil else {}))
    ax.axhline(TOL, color=ink2, lw=1, ls=(0, (4, 3)))
    ax.set_yscale("log"); ax.set_xlabel("Iterasi"); ax.set_ylabel("Selisih terbesar antar-iterasi (°C)")
    ax.grid(axis="y", color=grid, lw=0.8); ax.set_axisbelow(True)
    if kecil:
        ax.text(1, 1.6e-6, "batas berhenti 10⁻⁶", ha="left", va="bottom", color=ink2, fontsize=8.5)
        ax.text(len(hg) + 1, hg[-1] * 3.5, f"berhenti di\niterasi {len(hg)}", color=c1, fontsize=9, fontweight="bold", ha="left", va="bottom")
        ax.text(len(hj) + 1.2, 2.5e-6, f"berhenti di\niterasi {len(hj)}", color=c2, fontsize=9, fontweight="bold", ha="left", va="bottom")
        ax.set_xlim(0, len(hj) + 11)
    else:
        ax.text(60, 1.5e-6, "batas berhenti 10⁻⁶", ha="left", va="bottom", color=ink2, fontsize=8.5)
        ax.text(len(hg) + 150, 2.2e-6, f"Gauss-Seidel\n{len(hg):,} iterasi".replace(",", "."), color=c1, ha="left", va="bottom", fontsize=9, fontweight="bold")
        ax.text(len(hj) - 60, 2.2e-5, f"Jacobi\n{len(hj):,} iterasi".replace(",", "."), color=c2, ha="right", va="bottom", fontsize=9, fontweight="bold")
    ax.legend(frameon=False, loc="upper right")
    ax.set_title(f"Konvergensi contoh {n}×{n} (tebakan awal 0 °C)" if kecil else f"Konvergensi pada grid {n}×{n} (tebakan awal 0 °C)")
    simpan(fig, nama_file); plt.close(fig)


def gambar_peta_39x39():
    plt, cmap, ink, ink2, grid = siapkan_matplotlib()
    n = 39
    h = ambil(n)
    pad = padded_hasil(h["eliminasi"]["T"], n)
    fig, ax = plt.subplots(figsize=(5.4, 4.6), dpi=200)
    im = ax.imshow(pad, cmap=cmap, vmin=25, vmax=85, interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(n / 2 + 0.5, -2.2, "sisi atas 85 °C", ha="center", color=ink2, fontsize=9)
    ax.text(n / 2 + 0.5, n + 3.8, "sisi bawah 25 °C", ha="center", color=ink2, fontsize=9)
    ax.text(-2.2, n / 2 + 0.5, "sisi kiri 40 °C", rotation=90, va="center", ha="right", color=ink2, fontsize=9)
    ax.text(n + 3.2, n / 2 + 0.5, "sisi kanan 40 °C", rotation=270, va="center", ha="left", color=ink2, fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.1); cb.set_label("Suhu stabil (°C)"); cb.outline.set_visible(False)
    ax.set_title("Suhu stabil pada grid 39×39 (1.521 titik dalam)", pad=34)
    simpan(fig, "gambar3_peta_suhu_39x39.png"); plt.close(fig)


def menu_grafik():
    pilihan = {
        "1": ("Ilustrasi grid 3x3 (pelat, suhu tepi, T1-T9)", [gambar_pelat_3x3]),
        "2": ("Grafik konvergensi grid 3x3", [lambda: grafik_konvergensi(3, "gambar2_konvergensi_3x3.png")]),
        "3": ("Ilustrasi suhu grid 39x39", [gambar_peta_39x39]),
        "4": ("Grafik konvergensi grid 39x39", [lambda: grafik_konvergensi(39, "gambar4_konvergensi_39x39.png")]),
    }
    semua = [f for _, fs in pilihan.values() for f in fs]
    garis("=")
    print("GRAFIK / ILUSTRASI  (disimpan sebagai PNG di folder 'gambar_program')")
    garis("=")
    for k, (teks, _) in pilihan.items():
        print(f"  {k}. {teks}")
    print("  5. Semua gambar di atas")
    print("  0. Kembali")
    s = input("Pilih: ").strip()
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print("  matplotlib belum terpasang (pip install matplotlib).")
        return
    if s in pilihan:
        for f in pilihan[s][1]:
            f()
    elif s == "5":
        for f in semua:
            f()
    elif s != "0":
        print("  Pilihan tidak dikenal.")



# MENU 4: UJI KEBENARAN

def tiga_metode(n, bnd, init):
    A, b = build_system(n, bnd)
    return (gauss_elimination(A, b), gauss_seidel(n, bnd, init)[0].flatten(), jacobi(n, bnd, init)[0].flatten())


def terburuk(n, bnd, benar):
    w = [0.0, 0.0, 0.0]
    for init in (0, 25):
        for i, T in enumerate(tiga_metode(n, bnd, init)):
            w[i] = max(w[i], float(np.max(np.abs(T - benar))))
    return w


def menu_uji():
    garis("=")
    print("UJI KEBENARAN (error maksimum terhadap jawaban yang diketahui, terburuk dari tebakan awal 0 dan 25)")
    garis("=")
    xs = np.array([(c + 1) / 4 for c in range(3)])
    kasus = [
        ("1. Semua tepi 100 C (3x3) -> semua titik 100 C", 3, boundary(3, 100, 100, 100, 100), np.full(9, 100.0)),
        ("2. Kiri 0, kanan 100, atas-bawah naik linear (3x3) -> T = 100x", 3,
         boundary(3, lambda x: 100 * x, 0.0, 100.0, lambda x: 100 * x), np.tile(100 * xs, 3)),
        ("3. Contoh 2x2: atas 100, kiri 75, kanan 50, bawah 0", 2, boundary(2, 100, 75, 50, 0),
         np.array([71.875, 65.625, 46.875, 40.625])),
    ]
    for teks, n, bnd, benar in kasus:
        e, g, j = terburuk(n, bnd, benar)
        print(teks)
        print(f"     Eliminasi Gauss : {sci(e)}")
        print(f"     Gauss-Seidel    : {sci(g)}")
        print(f"     Jacobi          : {sci(j)}")
    garis("-", 82)
    A, b = build_system(3, boundary(3))
    print("\nVektor b 3x3 (125 85 125 40 0 40 65 25 65):", "sama dengan proposal" if list(b) == [125, 85, 125, 40, 0, 40, 65, 25, 65] else "BEDA")
    input("\nTekan Enter untuk kembali ke menu...")



# MENU UTAMA

def main():
    while True:
        print()
        garis("=")
        print("  DISTRIBUSI PANAS PADA CHIP  -  SPL (Gauss-Seidel, Jacobi, Eliminasi Gauss)")
        print("  Kelas A - Kelompok 3 - Metode Numerik")
        garis("=")
        print("  1. Grid 3x3    (9 persamaan)")
        print("  2. Grid 39x39  (1.521 persamaan)")
        print("  3. Grafik / ilustrasi")
        print("  4. Uji kebenaran")
        print("  0. Keluar")
        try:
            s = input("Pilih menu: ").strip()
        except EOFError:
            break
        try:
            if s == "1":
                menu_3x3()
            elif s == "2":
                menu_39x39()
            elif s == "3":
                menu_grafik()
            elif s == "4":
                menu_uji()
            elif s == "0":
                break
            else:
                print("Pilihan tidak dikenal.")
        except EOFError:
            break
    print("Selesai.")


if __name__ == "__main__":
    main()
