"""
Grafik dan gambar: ilustrasi grid dan grafik konvergensi (disimpan sebagai PNG).
matplotlib hanya dibutuhkan di berkas ini.
"""
import os

import numpy as np

from konfigurasi import TOP, SIDE, BOTTOM, TOL, FOLDER_GAMBAR
from penyelesaian import hitung_grid
from tabel import fx


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
    print(f"  tersimpan: {path}")


def grid_dengan_tepi(T, n):
    """Grid (n+2) x (n+2): titik dalam dari hasil T, sisi tepi diisi suhu sisi."""
    pad = np.full((n + 2, n + 2), np.nan)
    pad[1:-1, 1:-1] = np.asarray(T).reshape(n, n)
    pad[0, 1:-1], pad[-1, 1:-1], pad[1:-1, 0], pad[1:-1, -1] = TOP, BOTTOM, SIDE, SIDE
    return pad


def gambar_pelat_3x3():
    plt, cmap, ink, ink2, grid = siapkan_matplotlib()
    h = hitung_grid(3)
    pad = grid_dengan_tepi(h["eliminasi"]["T"], 3)
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
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title("Suhu stabil contoh 3×3 (sel tepi = suhu sisi yang diketahui)", pad=10)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.04)
    cb.set_label("Suhu (°C)")
    cb.outline.set_visible(False)
    simpan(fig, "gambar1_pelat_3x3.png")
    plt.close(fig)


def grafik_konvergensi(n, nama_file):
    """Selisih maks per iterasi sampai konvergen (grid 3x3 atau 39x39), skala log."""
    plt, cmap, ink, ink2, grid = siapkan_matplotlib()
    c1, c2 = "#2a78d6", "#eb6834"      # biru = Gauss-Seidel, oranye = Jacobi
    h = hitung_grid(n)
    hg, hj = h["gs"]["riwayat"], h["jacobi"]["riwayat"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6), dpi=200)
    kecil = n <= 3
    ax.plot(range(1, len(hj) + 1), hj, color=c2, lw=2, label="Jacobi", **({"marker": "o", "ms": 3.5} if kecil else {}))
    ax.plot(range(1, len(hg) + 1), hg, color=c1, lw=2, label="Gauss-Seidel", **({"marker": "o", "ms": 3.5} if kecil else {}))
    ax.axhline(TOL, color=ink2, lw=1, ls=(0, (4, 3)))
    ax.set_yscale("log")
    ax.set_xlabel("Iterasi")
    ax.set_ylabel("Selisih terbesar antar-iterasi (°C)")
    ax.grid(axis="y", color=grid, lw=0.8)
    ax.set_axisbelow(True)
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
    simpan(fig, nama_file)
    plt.close(fig)


def gambar_peta_39x39():
    plt, cmap, ink, ink2, grid = siapkan_matplotlib()
    n = 39
    h = hitung_grid(n, singkat=False)
    pad = grid_dengan_tepi(h["eliminasi"]["T"], n)
    fig, ax = plt.subplots(figsize=(5.4, 4.6), dpi=200)
    im = ax.imshow(pad, cmap=cmap, vmin=25, vmax=85, interpolation="nearest")
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.text(n / 2 + 0.5, -2.2, "sisi atas 85 °C", ha="center", color=ink2, fontsize=9)
    ax.text(n / 2 + 0.5, n + 3.8, "sisi bawah 25 °C", ha="center", color=ink2, fontsize=9)
    ax.text(-2.2, n / 2 + 0.5, "sisi kiri 40 °C", rotation=90, va="center", ha="right", color=ink2, fontsize=9)
    ax.text(n + 3.2, n / 2 + 0.5, "sisi kanan 40 °C", rotation=270, va="center", ha="left", color=ink2, fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.1)
    cb.set_label("Suhu stabil (°C)")
    cb.outline.set_visible(False)
    ax.set_title("Suhu stabil pada grid 39×39 (1.521 titik dalam)", pad=34)
    simpan(fig, "gambar3_peta_suhu_39x39.png")
    plt.close(fig)


PILIHAN = {
    "1": ("Ilustrasi grid 3x3 (pelat, suhu tepi, T1-T9)", [gambar_pelat_3x3]),
    "2": ("Grafik konvergensi grid 3x3", [lambda: grafik_konvergensi(3, "gambar2_konvergensi_3x3.png")]),
    "3": ("Ilustrasi suhu grid 39x39", [gambar_peta_39x39]),
    "4": ("Grafik konvergensi grid 39x39 (sampai konvergen)", [lambda: grafik_konvergensi(39, "gambar4_konvergensi_39x39.png")]),
}


def buat(pilih):
    """Buat gambar sesuai pilihan '1'..'4' atau '5' (semua)."""
    if pilih == "5":
        daftar = [f for _, fs in PILIHAN.values() for f in fs]
    elif pilih in PILIHAN:
        daftar = PILIHAN[pilih][1]
    else:
        print("  Pilihan tidak dikenal.")
        return
    for f in daftar:
        f()
