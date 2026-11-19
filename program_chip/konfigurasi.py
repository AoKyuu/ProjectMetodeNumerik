"""Konfigurasi: konstanta kasus, syarat berhenti, dan folder keluaran."""
import os

FOLDER_PROGRAM = os.path.dirname(os.path.abspath(__file__))
FOLDER_GAMBAR = os.path.join(FOLDER_PROGRAM, "gambar_program")   # hasil menu 3
FOLDER_HASIL = os.path.join(FOLDER_PROGRAM, "hasil_program")     # hasil simpan CSV

# Suhu sisi pelat (derajat Celsius)
TOP, SIDE, BOTTOM = 85.0, 40.0, 25.0     # atas, kiri = kanan, bawah

# Metode iteratif
TOL = 1e-6              # berhenti jika selisih terbesar antar-iterasi < TOL
MAXIT = 200000          # batas pengaman jumlah iterasi
ITER_SINGKAT = 50       # jumlah iterasi pada mode singkat (grid 39x39)
