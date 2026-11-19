"""
Tabel dan hasil: semua yang ditampilkan ke layar (tabel, ringkasan, langkah perhitungan) dan penyimpanan CSV.
Angka ditulis dengan koma desimal dan tanpa notasi ilmiah.
"""
import os
from decimal import Decimal

import numpy as np

from konfigurasi import FOLDER_HASIL
from penyelesaian import tetangga


# ------------------------------------------------------------------ format angka
def fx(x, d=4):
    """Angka dengan d angka di belakang koma, memakai koma desimal."""
    return f"{x:.{d}f}".replace(".", ",")


def angka_kecil(x):
    """Angka kecil ditulis biasa (tanpa huruf e), 2 angka penting, mis. 0,00000059; angka besar 4 desimal."""
    if abs(x) >= 0.01:
        return fx(x, 4)
    return format(Decimal(f"{x:.1e}"), "f").replace(".", ",")


def waktu_teks(detik):
    if detik < 1e-3:
        return f"{detik * 1e6:,.0f} mikrodetik".replace(",", ".")
    return f"{detik:.3f} detik".replace(".", ",")


def konstanta(v):
    return str(int(v)) if float(v) == int(v) else fx(v, 7)


def garis(ch="-", n=78):
    print(ch * n)


def judul(teks):
    garis("=")
    print(teks)
    garis("=")


# ------------------------------------------------------------------ sistem persamaan 3x3
def tampil_sistem(A, b):
    N = len(b)
    nb = tetangga(A)
    print("\nPersamaan tiap titik dalam (T = rata-rata 4 tetangga, suhu tepi sudah dimasukkan):")
    for i in range(N):
        suku = ([konstanta(b[i])] if b[i] != 0 else []) + [f"T{j + 1}" for j in nb[i]]
        print(f"  T{i + 1} = ({' + '.join(suku)}) / 4")
    print(f"\nMatriks A ({N}x{N}):")
    for i in range(N):
        print("  [" + " ".join(f"{int(v):>2}" for v in A[i]) + " ]")
    print("Vektor b:", [int(v) for v in b])


def tampil_hasil(h):
    """Tabel suhu T1..TN dari tiga metode (hanya untuk grid kecil)."""
    N = h["N"]
    print("\nHasil suhu stabil (derajat C):")
    kepala = f"{'Titik':<7}{'Gauss-Seidel':>16}{'Jacobi':>16}{'Eliminasi Gauss':>18}"
    print(kepala)
    garis("-", len(kepala))
    for i in range(N):
        print(f"T{i + 1:<6}{fx(h['gs']['T'][i]):>16}{fx(h['jacobi']['T'][i]):>16}{fx(h['eliminasi']['T'][i]):>18}")
    garis("-", len(kepala))


def ringkasan(h):
    """Ringkasan jumlah iterasi, waktu, error, dan residual."""
    print(f"\n{'Metode':<17}{'Iterasi':>8}{'Waktu':>19}{'Error maks (C)':>19}{'Residual maks':>19}")
    garis("-", 82)
    for nama, kunci in (("Gauss-Seidel", "gs"), ("Jacobi", "jacobi")):
        r = h[kunci]
        print(f"{nama:<17}{r['iterasi']:>8}{waktu_teks(r['waktu']):>19}{angka_kecil(r['error']):>19}{angka_kecil(r['residual']):>19}")
    e = h["eliminasi"]
    print(f"{'Eliminasi Gauss':<17}{'langsung':>8}{waktu_teks(e['waktu']):>19}{'(acuan)':>19}{angka_kecil(e['residual']):>19}")
    garis("-", 82)
    print("Error = selisih mutlak terbesar terhadap Eliminasi Gauss.  Residual = nilai mutlak terbesar A.T - b.")
    if h["singkat"]:
        print(f"Gauss-Seidel dan Jacobi dihentikan pada iterasi {h['gs']['iterasi']} (mode singkat), belum konvergen.")
    else:
        print("Syarat berhenti: selisih terbesar antar-iterasi < 0,000001 (bukan jaminan error < 0,000001).")


# ------------------------------------------------------------------ metode iteratif: langkah dan tabel (grid kecil)
NAMA_METODE = {"gs": "Gauss-Seidel", "jacobi": "Jacobi"}


def rumus_iterasi(A, b, metode):
    """Cetak rumus tiap titik dengan penanda (k) dan (k-1)."""
    N = len(b)
    nb = tetangga(A)
    print(f"\nRumus iterasi {NAMA_METODE[metode]}:")
    for i in range(N):
        suku = [konstanta(b[i])] if b[i] != 0 else []
        for j in nb[i]:
            baru = (metode == "gs" and j < i)
            suku.append(f"T{j + 1}({'k' if baru else 'k-1'})")
        print(f"  T{i + 1}(k) = ({' + '.join(suku)}) / 4")


def langkah_iterasi(h, metode, jumlah=5):
    """Perhitungan iterasi 1..jumlah: tiap titik ditulis dengan angka (7 desimal)."""
    A, b, r = h["A"], h["b"], h[metode]
    N = len(b)
    nb = tetangga(A)
    rumus_iterasi(A, b, metode)
    print(f"\nPerhitungan iterasi 1 sampai {jumlah} (tebakan awal {fx(h['init'], 1)}), 7 angka di belakang koma:")
    for k in range(1, min(jumlah, r["iterasi"]) + 1):
        cur = r["snaps"][k - 1]
        lama = r["snaps"][k - 2] if k > 1 else np.full(N, h["init"])
        print(f"\nIterasi {k}" + (f" (k = 1, semua nilai iterasi k-1 = {fx(h['init'], 0)})" if k == 1 else ""))
        for i in range(N):
            nilai = [cur[j] if (metode == "gs" and j < i) else lama[j] for j in nb[i]]
            suku = ([konstanta(b[i])] if b[i] != 0 else []) + [fx(v, 7) for v in nilai]
            print(f"  T{i + 1} = ({' + '.join(suku)}) / 4 = {fx(cur[i], 7)}")
        print(f"  Selisih maks = {fx(r['riwayat'][k - 1], 7)}")
    print("\nCatatan: angka dihitung dengan presisi penuh lalu ditampilkan 7 desimal; bila dihitung ulang dari angka\n"
          "yang sudah dibulatkan, beberapa sel dapat bergeser 0,0000001.")


def tabel_iterasi(h, metode):
    """Tabel semua iterasi (0 sampai berhenti) dengan 7 angka di belakang koma."""
    r = h[metode]
    N = h["N"]
    print(f"\nTabel iterasi {NAMA_METODE[metode]} (tebakan awal {fx(h['init'], 1)} C), 7 angka di belakang koma")
    kepala = "Iter " + "".join(f"{'T' + str(i + 1):>13}" for i in range(N)) + f"{'Selisih maks':>14}"
    print(kepala)
    garis("-", len(kepala))
    print(f"{0:<5}" + "".join(f"{fx(h['init'], 7):>13}" for _ in range(N)))
    for i, (snap, d) in enumerate(zip(r["snaps"], r["riwayat"]), start=1):
        print(f"{i:<5}" + "".join(f"{fx(v, 7):>13}" for v in snap) + f"{fx(d, 7):>14}")
    print(f"Berhenti pada iterasi {r['iterasi']}.")


# ------------------------------------------------------------------ eliminasi Gauss: langkah demi langkah (grid kecil)
def cetak_matriks(M, judul_matriks, ada_c=True, semua_desimal=False):
    N = M.shape[0]
    print(f"\n{judul_matriks}")
    kepala = f"{'':<5}" + "".join(f"{'T' + str(j + 1):>12}" for j in range(N)) + (f"{'c':>13}" if ada_c else "")
    print(kepala)
    garis("-", len(kepala))
    for i in range(N):
        sel = []
        for j in range(M.shape[1]):
            v = M[i, j]
            if abs(v) < 1e-12:
                t = "0"
            elif not semua_desimal and abs(v - round(v)) < 1e-12:
                t = str(int(round(v)))
            else:
                t = fx(v, 7)
            sel.append(f"{t:>12}" if j < N else f"{t:>13}")
        print(f"R{i + 1:<4}" + "".join(sel))


def langkah_eliminasi(info, tampil_tiap_langkah=False):
    """Cetak [A|b], operasi baris tiap langkah, [U|c], dan substitusi mundur (rumus lalu angka)."""
    U, c, x = info["U"], info["c"], info["x"]
    N = len(c)
    cetak_matriks(info["awal"], "Matriks augmented [A | b]")
    print("\nEliminasi maju: m = a(i,k) / a(k,k);  R(i) <- R(i) - m * R(k)")
    for st in info["langkah"]:
        k = st["k"]
        if st["tukar"]:
            print(f"\nLangkah {k + 1}: tukar R{st['tukar'][0] + 1} dengan R{st['tukar'][1] + 1} (pivoting)")
        else:
            print(f"\nLangkah {k + 1}")
        print(f"  pivot a{k + 1}{k + 1} = {fx(st['pivot'], 7)}")
        for (i, aik, m) in st["ops"]:
            print(f"  m{i + 1}{k + 1} = a{i + 1}{k + 1} / a{k + 1}{k + 1} = ({fx(aik, 7)}) / {fx(st['pivot'], 7)} = {fx(m, 7)};"
                  f"   R{i + 1} <- R{i + 1} - ({fx(m, 7)}) * R{k + 1}")
        if tampil_tiap_langkah:
            cetak_matriks(st["M"], f"Matriks setelah langkah {k + 1}")
    cetak_matriks(np.hstack([U, c[:, None]]), "Matriks segitiga atas [U | c]", semua_desimal=True)
    print("\nSubstitusi mundur (U.T = c, dari baris terakhir ke baris 1):")
    for i in range(N - 1, -1, -1):
        js = [j for j in range(i + 1, N) if abs(U[i, j]) > 1e-12]
        if js:
            simb = f"T{i + 1} = (c{i + 1} - " + " - ".join(f"u{i + 1}{j + 1}*T{j + 1}" for j in js) + f") / u{i + 1}{i + 1}"
            angka = (f"T{i + 1} = ({fx(c[i], 7)} - " + " - ".join(f"({fx(U[i, j], 7)})({fx(x[j], 7)})" for j in js)
                     + f") / {fx(U[i, i], 7)} = {fx(x[i], 7)}")
        else:
            simb = f"T{i + 1} = c{i + 1} / u{i + 1}{i + 1}"
            angka = f"T{i + 1} = {fx(c[i], 7)} / {fx(U[i, i], 7)} = {fx(x[i], 7)}"
        print(f"  {simb}")
        print(f"  {angka}")
    print("\nHasil Eliminasi Gauss (7 angka di belakang koma):")
    for i in range(N):
        print(f"  T{i + 1} = {fx(x[i], 7)}")


def hasil_eliminasi_saja(h):
    """Hanya hasil akhir Eliminasi Gauss (tanpa langkah)."""
    e = h["eliminasi"]
    print(f"\nHasil akhir Eliminasi Gauss (waktu {waktu_teks(e['waktu'])}, residual maks {angka_kecil(e['residual'])}):")
    for i, v in enumerate(e["T"]):
        print(f"  T{i + 1} = {fx(v, 7)}")


# ------------------------------------------------------------------ grid besar: iterasi singkat dan cuplikan hasil
def tabel_singkat(h):
    """Ringkasan tiap iterasi untuk grid besar: selisih maks dan suhu titik tengah (tanpa menampilkan semua titik)."""
    n = h["n"]
    tengah = (n // 2) * n + (n // 2)
    gs, jc = h["gs"], h["jacobi"]
    print(f"\nIterasi singkat grid {n}x{n} (tebakan awal {fx(h['init'], 1)} C): selisih maks dan suhu titik tengah "
          f"(baris {n // 2 + 1}, kolom {n // 2 + 1})")
    kepala = f"{'Iter':<6}{'GS selisih maks':>17}{'GS T tengah':>14}{'Jacobi selisih maks':>21}{'Jacobi T tengah':>17}"
    print(kepala)
    garis("-", len(kepala))
    for k in range(1, gs["iterasi"] + 1):
        print(f"{k:<6}{fx(gs['riwayat'][k - 1], 7):>17}{fx(gs['snaps'][k - 1][tengah], 7):>14}"
              f"{fx(jc['riwayat'][k - 1], 7):>21}{fx(jc['snaps'][k - 1][tengah], 7):>17}")
    garis("-", len(kepala))


def cuplikan(T, n, d=2):
    """Cuplikan hasil akhir: suhu pada 7 baris x 7 kolom yang tersebar merata (termasuk titik tengah)."""
    idx = [0, n // 6, 2 * (n // 6), n // 2, n - 1 - 2 * (n // 6), n - 1 - (n // 6), n - 1]
    G = np.asarray(T).reshape(n, n)
    print(f"\nCuplikan hasil akhir Eliminasi Gauss ({len(idx)}x{len(idx)} titik dari {n}x{n}), suhu dalam C")
    label = "baris\\kolom"
    kepala = f"{label:<12}" + "".join(f"{j + 1:>9}" for j in idx)
    print(kepala)
    garis("-", len(kepala))
    for i in idx:
        print(f"{i + 1:<12}" + "".join(f"{fx(G[i, j], d):>9}" for j in idx))
    garis("-", len(kepala))
    print(f"Suhu terendah {fx(G.min(), 4)} C, tertinggi {fx(G.max(), 4)} C "
          f"(di dalam rentang 25-85 C: {'ya' if 25 <= G.min() and G.max() <= 85 else 'tidak'}); "
          f"titik tengah {fx(G[n // 2, n // 2], 4)} C")


def simpan_csv(T, n, nama="hasil_grid_39x39.csv"):
    """Simpan seluruh suhu hasil akhir (n x n) ke CSV: pemisah titik koma, koma desimal (cocok untuk Excel Indonesia)."""
    os.makedirs(FOLDER_HASIL, exist_ok=True)
    path = os.path.join(FOLDER_HASIL, nama)
    G = np.asarray(T).reshape(n, n)
    with open(path, "w", encoding="utf-8-sig") as f:
        f.write("baris\\kolom;" + ";".join(str(j + 1) for j in range(n)) + "\n")
        for i in range(n):
            f.write(f"{i + 1};" + ";".join(fx(G[i, j], 7) for j in range(n)) + "\n")
    print(f"  tersimpan: {path}")
    return path
