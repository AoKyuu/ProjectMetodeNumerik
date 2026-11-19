"""
Distribusi Panas pada Chip - Sistem Persamaan Linear (SPL)
Kelas A - Kelompok 3 - Metode Numerik (Tahap 2)

Jalankan:  python main.py

Susunan berkas:
  main.py          antarmuka dan pilihan menu di terminal (berkas ini)
  penyelesaian.py  model A.T = b dan tiga metode (Gauss-Seidel, Jacobi, Eliminasi Gauss)
  tabel.py         tabel, ringkasan, langkah perhitungan, dan penyimpanan CSV
  grafik.py        ilustrasi dan grafik (PNG), membutuhkan matplotlib
  uji.py           uji kebenaran (kasus dengan jawaban diketahui)
  konfigurasi.py   suhu sisi, syarat berhenti, dan folder keluaran
"""
import penyelesaian as pe
import tabel as tb
import grafik
import uji


def tanya_tebakan():
    s = input("Tebakan awal semua titik dalam (derajat C) [Enter = 0]: ").strip().replace(",", ".")
    if s == "":
        return 0.0
    try:
        return float(s)
    except ValueError:
        print("  Input bukan angka, dipakai 0.")
        return 0.0


def tunggu():
    input("\nTekan Enter untuk melanjutkan...")


# =====================================================================================
# MENU 1: GRID 3x3
# =====================================================================================
def menu_3x3():
    tb.judul("GRID 3x3  (9 titik dalam, 9 persamaan, matriks A 9x9)")
    print("Suhu sisi: atas 85 C, kiri dan kanan 40 C, bawah 25 C")
    init = tanya_tebakan()
    h = pe.hitung_grid(3, init)
    while True:
        print("\n--- Menu grid 3x3 ---")
        print("  1. Ringkasan tiga metode (persamaan, matriks, hasil T1-T9, perbandingan)")
        print("  2. Gauss-Seidel : iterasi 1-5 per langkah + tabel semua iterasi")
        print("  3. Jacobi       : iterasi 1-5 per langkah + tabel semua iterasi")
        print("  4. Eliminasi Gauss : per langkah (matriks augmented, operasi baris, substitusi mundur)")
        print("  5. Eliminasi Gauss : hasil akhir saja")
        print("  0. Kembali")
        s = input("Pilih: ").strip()
        if s == "1":
            tb.tampil_sistem(h["A"], h["b"])
            tb.tampil_hasil(h)
            tb.ringkasan(h)
            tunggu()
        elif s in ("2", "3"):
            metode = "gs" if s == "2" else "jacobi"
            tb.langkah_iterasi(h, metode, jumlah=5)
            tb.tabel_iterasi(h, metode)
            tunggu()
        elif s == "4":
            info = pe.eliminasi_langkah(h["A"], h["b"])
            tiap = input("Tampilkan matriks setelah tiap langkah? (y/N): ").strip().lower() == "y"
            tb.langkah_eliminasi(info, tampil_tiap_langkah=tiap)
            tunggu()
        elif s == "5":
            tb.hasil_eliminasi_saja(h)
            tunggu()
        elif s == "0" or s == "":
            return
        else:
            print("  Pilihan tidak dikenal.")


# =====================================================================================
# MENU 2: GRID 39x39
# =====================================================================================
def menu_39x39():
    tb.judul("GRID 39x39  (1.521 titik dalam, 1.521 persamaan, matriks A 1.521 x 1.521)")
    print("Suhu sisi: atas 85 C, kiri dan kanan 40 C, bawah 25 C")
    print("Suhu T1 sampai T1521 tidak ditabelkan (terlalu banyak); yang ditampilkan ringkasan dan cuplikan.")
    init = tanya_tebakan()
    while True:
        print("\n--- Menu grid 39x39 ---")
        print("  1. Tiga metode: iterasi singkat (50 iterasi) + Eliminasi Gauss langsung hasil akhir")
        print("  2. Gauss-Seidel dan Jacobi sampai konvergen (ringkasan) + Eliminasi Gauss")
        print("  3. Simpan hasil akhir 39x39 (Eliminasi Gauss) ke berkas CSV")
        print("  0. Kembali")
        s = input("Pilih: ").strip()
        if s == "1":
            h = pe.hitung_grid(39, init, singkat=True)
            tb.tabel_singkat(h)
            tb.ringkasan(h)
            tb.cuplikan(h["eliminasi"]["T"], 39)
            print("\nSetelah 50 iterasi kedua metode iteratif belum konvergen (lihat kolom Error maks).")
            print("Untuk sampai konvergen pilih menu 2.")
            tunggu()
        elif s == "2":
            h = pe.hitung_grid(39, init, singkat=False)
            tb.ringkasan(h)
            tb.cuplikan(h["eliminasi"]["T"], 39)
            tunggu()
        elif s == "3":
            h = pe.hitung_grid(39, init, singkat=True)
            tb.simpan_csv(h["eliminasi"]["T"], 39)
            tunggu()
        elif s == "0" or s == "":
            return
        else:
            print("  Pilihan tidak dikenal.")


# =====================================================================================
# MENU 3: GRAFIK / ILUSTRASI
# =====================================================================================
def menu_grafik():
    tb.judul("GRAFIK / ILUSTRASI  (disimpan sebagai PNG di folder 'gambar_program')")
    for k, (teks, _) in grafik.PILIHAN.items():
        print(f"  {k}. {teks}")
    print("  5. Semua gambar di atas")
    print("  0. Kembali")
    s = input("Pilih: ").strip()
    if s in ("", "0"):
        return
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        print("  matplotlib belum terpasang (pip install matplotlib).")
        return
    grafik.buat(s)
    tunggu()


# =====================================================================================
# MENU 4: UJI KEBENARAN
# =====================================================================================
def menu_uji():
    tb.judul("UJI KEBENARAN (error maksimum terhadap jawaban yang diketahui, terburuk dari tebakan awal 0 dan 25)")
    uji.jalankan_uji()
    tunggu()


# =====================================================================================
# MENU UTAMA
# =====================================================================================
def main():
    while True:
        print()
        tb.garis("=")
        print("  DISTRIBUSI PANAS PADA CHIP  -  SPL (Gauss-Seidel, Jacobi, Eliminasi Gauss)")
        print("  Kelas A - Kelompok 3 - Metode Numerik")
        tb.garis("=")
        print("  1. Grid 3x3    (9 persamaan)")
        print("  2. Grid 39x39  (1.521 persamaan)")
        print("  3. Grafik / ilustrasi")
        print("  4. Uji kebenaran")
        print("  0. Keluar")
        try:
            s = input("Pilih menu: ").strip()
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
        except (EOFError, KeyboardInterrupt):
            break
    print("Selesai.")


if __name__ == "__main__":
    main()
