"""Uji kebenaran: tiga kasus yang jawabannya diketahui, dijalankan pada tiga metode."""
import numpy as np

from penyelesaian import boundary, build_system, gauss_elimination, gauss_seidel, jacobi
from tabel import angka_kecil, garis

B_PROPOSAL = [125, 85, 125, 40, 0, 40, 65, 25, 65]


def tiga_metode(n, bnd, init):
    A, b = build_system(n, bnd)
    return (gauss_elimination(A, b), gauss_seidel(n, bnd, init)[0].flatten(), jacobi(n, bnd, init)[0].flatten())


def terburuk(n, bnd, benar):
    """Error maksimum tiap metode, diambil yang terburuk dari tebakan awal 0 dan 25."""
    w = [0.0, 0.0, 0.0]
    for init in (0, 25):
        for i, T in enumerate(tiga_metode(n, bnd, init)):
            w[i] = max(w[i], float(np.max(np.abs(T - benar))))
    return w


def jalankan_uji():
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
        print(f"     Eliminasi Gauss : {angka_kecil(e)}")
        print(f"     Gauss-Seidel    : {angka_kecil(g)}")
        print(f"     Jacobi          : {angka_kecil(j)}")
    garis("-", 82)
    _, b = build_system(3, boundary(3))
    print("Vektor b 3x3 (125 85 125 40 0 40 65 25 65):", "sama dengan proposal" if list(b) == B_PROPOSAL else "BEDA")
