# Simulasi Distribusi Panas pada Chip

**Laporan Tugas Akhir Mata Kuliah Metode Numerik**  
**Kelas A — Kelompok 3**  
Program Studi Informatika PSDKU Kebumen  
Universitas Sebelas Maret — 2026

## 1\. Deskripsi Program

Program Python ini menyelesaikan model distribusi suhu stabil (*steady-state*) pada bagian dalam chip prosesor yang dimodelkan sebagai pelat persegi tipis. Sistem disusun sebagai Sistem Persamaan Linear (SPL) dan diselesaikan menggunakan tiga metode numerik:

* **Gauss-Seidel**, sebagai metode iteratif utama.
* **Jacobi**, sebagai pembanding metode iteratif.
* **Eliminasi Gauss dengan pivoting parsial**, sebagai metode langsung dan solusi acuan.

Kasus dan parameter mengikuti proposal Tahap 1 dan laporan tugas akhir. Program menyediakan simulasi grid 3×3 dan 39×39, grafik/ilustrasi, serta menu uji kebenaran.

## 2\. Model dan Parameter

Pada keadaan stabil, suhu setiap titik dalam sama dengan rata-rata suhu empat tetangganya: atas, bawah, kiri, dan kanan. Untuk titik dalam pada baris (i) dan kolom (j):

<img width="457" height="94" alt="image" src="https://github.com/user-attachments/assets/ca710332-32d6-4cca-9032-d8b2e1851447" />


Suhu pada sisi pelat yang diketahui dimasukkan ke ruas kanan. Persamaan seluruh titik dalam kemudian disusun menjadi:

<img width="135" height="58" alt="image" src="https://github.com/user-attachments/assets/87d72884-e6a7-4f63-9a4f-34157b8bffa1" />


Parameter suhu batas yang digunakan:

|Sisi pelat|Suhu batas|
|-|-:|
|Atas|85 °C|
|Kiri|40 °C|
|Kanan|40 °C|
|Bawah|25 °C|

Setiap titik dalam menghasilkan satu persamaan. Penomoran titik dilakukan per baris dari kiri atas, yaitu T1 sampai T9 untuk grid 3×3. Matriks grid 3×3 berukuran 9×9, sedangkan matriks grid 39×39 berukuran 1.521×1.521. Matriks memiliki paling banyak lima elemen tidak nol pada setiap baris.

Parameter komputasi:

|Parameter|Nilai|Keterangan|
|-|-:|-|
|`TOP`|`85.0`|Suhu sisi atas (°C)|
|`SIDE`|`40.0`|Suhu sisi kiri dan kanan (°C)|
|`BOTTOM`|`25.0`|Suhu sisi bawah (°C)|
|`TOL`|`1e-6`|Batas selisih maksimum antariterasi|
|`MAXIT`|`200000`|Batas maksimum iterasi|

## 3\. Metode Numerik

### Gauss-Seidel

Titik dalam diperbarui berurutan dari T1 sampai T9. Nilai tetangga yang sudah diperbarui pada iterasi yang sama langsung digunakan, sedangkan nilai tetangga yang belum diperbarui masih menggunakan nilai iterasi sebelumnya.

### Jacobi

Seluruh titik pada iterasi baru dihitung menggunakan nilai dari iterasi sebelumnya. Oleh karena itu, metode ini menggunakan dua salinan grid.

### Eliminasi Gauss

SPL diselesaikan melalui eliminasi maju dengan pivoting parsial, kemudian dilanjutkan dengan substitusi mundur. Solusinya digunakan sebagai acuan untuk menghitung error Gauss-Seidel dan Jacobi. Hasil juga diperiksa silang dengan `numpy.linalg.solve` selisih maksimum yang dilaporkan adalah sekitar (1,4\\times10^{-14}) °C.

### Pengaturan Percobaan

* Tebakan awal yang diuji: seluruh titik dalam bernilai 0 °C dan 25 °C.
* Syarat berhenti Gauss-Seidel dan Jacobi: selisih terbesar antariterasi pada seluruh titik kurang dari (10^{-6}).
* Batas maksimum iterasi: 200.000.
* Error dihitung sebagai selisih mutlak maksimum hasil metode iteratif terhadap solusi Eliminasi Gauss.
* Residual dihitung sebagai nilai mutlak maksimum dari (A\\mathbf{T}-\\mathbf{b}).

## 4\. Persyaratan dan Instalasi

Persyaratan:

* Python 3.
* NumPy.
* Matplotlib (digunakan untuk menu grafik/ilustrasi).

Pasang dependensi dari folder proyek:

```bash
python -m pip install -r NumPy

python -m pip install -r Matplotlib

```

Untuk menjalankan program, gunakan IDE yang mendukung bahasa pemrograman Python.



Saat menu grid meminta tebakan awal, tekan **Enter** untuk menggunakan 0 °C. Pengguna juga dapat memasukkan 25 untuk menguji tebakan awal 25 °C. Input desimal menerima titik maupun koma.

## 5\. Menu Program

|Menu|Fitur|Keluaran utama|
|-|-|-|
|`1`|Grid 3×3 (9 persamaan)|Meminta tebakan awal; menampilkan sembilan persamaan, matriks A 9×9, vektor b, hasil ketiga metode, ringkasan iterasi/waktu/error/residual, serta pilihan menampilkan tabel seluruh iterasi Gauss-Seidel atau Jacobi.|
|`2`|Grid 39×39 (1.521 persamaan)|Meminta tebakan awal; menjalankan ketiga metode dan menampilkan ringkasan hasil, suhu titik tengah, suhu minimum, dan suhu maksimum. Suhu T1 sampai T1521 tidak ditampilkan satu per satu.|
|`3`|Grafik / ilustrasi|Membuat empat gambar PNG: ilustrasi grid 3×3, grafik konvergensi 3×3, ilustrasi suhu grid 39×39, dan grafik konvergensi 39×39. Pengguna dapat memilih satu gambar atau semua gambar.|
|`4`|Uji kebenaran|Menjalankan tiga kasus yang jawabannya diketahui dengan tebakan awal 0 °C dan 25 °C, menampilkan error maksimum, serta memeriksa vektor b grid 3×3 terhadap proposal.|
|`0`|Keluar|Mengakhiri program.|

## 6\. Hasil yang Dilaporkan

### 6.1 Grid 3×3

Hasil akhir suhu untuk tebakan awal 0 °C ditampilkan berikut (°C):

|Titik|Gauss-Seidel (27 iterasi)|Jacobi (50 iterasi)|Eliminasi Gauss|
|-|-:|-:|-:|
|T1|58,2143|58,2143|58,2143|
|T2|62,2321|62,2321|62,2321|
|T3|58,2143|58,2143|58,2143|
|T4|45,6250|45,6250|45,6250|
|T5|47,5000|47,5000|47,5000|
|T6|45,6250|45,6250|45,6250|
|T7|36,7857|36,7857|36,7857|
|T8|36,5179|36,5179|36,5179|
|T9|36,7857|36,7857|36,7857|

Jumlah iterasi sampai selisih maksimum antariterasi kurang dari (10^{-6}):

|Metode|Tebakan awal 0 °C|Tebakan awal 25 °C|
|-|-:|-:|
|Gauss-Seidel|27|25|
|Jacobi|50|48|
|Rasio Jacobi / Gauss-Seidel|1,85|1,92|
|Eliminasi Gauss|Langsung (tanpa iterasi)|Langsung (tanpa iterasi)|

Gauss-Seidel membutuhkan sekitar separuh jumlah iterasi Jacobi pada kedua tebakan awal. Tebakan awal 25 °C mengurangi jumlah iterasi sebanyak dua untuk masing-masing metode iteratif dan tidak mengubah hasil akhir.

### 6.2 Waktu Komputasi Grid 3×3

Waktu berikut adalah waktu terbaik dari 500 kali eksekusi pada pengujian laporan. Hasil dapat berbeda menurut perangkat dan cara penulisan kode.

|Gauss-Seidel|Jacobi|Eliminasi Gauss|
|-:|-:|-:|
|31 mikrodetik|60 mikrodetik|49 mikrodetik|

Karena hanya terdapat sembilan persamaan, waktu ketiga metode sangat kecil dan tidak digunakan sebagai dasar utama kesimpulan perbandingan.

### 6.3 Akurasi Grid 3×3

Error maksimum terhadap solusi Eliminasi Gauss dan residual untuk tebakan awal 0 °C:

|Metode|Error maksimum (awal 0 °C)|Error maksimum (awal 25 °C)|Residual maksimum (awal 0 °C)|
|-|-:|-:|-:|
|Gauss-Seidel|5,9×10⁻⁷ °C|9,8×10⁻⁷ °C|1,2×10⁻⁶|
|Jacobi|2,1×10⁻⁶ °C|2,0×10⁻⁶ °C|2,8×10⁻⁶|
|Eliminasi Gauss|Solusi acuan|Solusi acuan|2,8×10⁻¹⁴|

Eliminasi Gauss menghasilkan residual yang sangat kecil dan digunakan sebagai solusi acuan. Error Jacobi sedikit lebih besar daripada batas perubahan antariterasi (10^{-6}), sesuai dengan sifat kriteria berhenti yang digunakan.

### 6.4 Uji Kebenaran

Ketiga kasus berikut diuji menggunakan tebakan awal 0 °C dan 25 °C. Tabel menampilkan error maksimum terburuk terhadap jawaban yang diketahui.

|Kasus|Jawaban yang diharapkan|Gauss-Seidel|Jacobi|Eliminasi Gauss|
|-|-|-:|-:|-:|
|Semua tepi 100 °C (grid 3×3)|Semua titik 100 °C|7,0×10⁻⁷|2,2×10⁻⁶|2,8×10⁻¹⁴|
|Kiri 0 °C, kanan 100 °C, atas-bawah naik linear (grid 3×3)|(T=100x)|8,5×10⁻⁷|2,2×10⁻⁶|1,4×10⁻¹⁴|
|Contoh batas atas 100 °C, kiri 75 °C, kanan 50 °C, bawah 0 °C (grid 2×2)|T1=71,875; T2=65,625; T3=46,875; T4=40,625|2,3×10⁻⁷|9,3×10⁻⁷|7,1×10⁻¹⁵|

Ketiga metode lolos pada seluruh kasus uji. Error terbesar yang dilaporkan adalah sekitar (2,2\\times10^{-6}) °C. Vektor b grid 3×3 pada konfigurasi standar juga sama dengan proposal:

```text
\[125, 85, 125, 40, 0, 40, 65, 25, 65]
```

Pemeriksaan solusi Eliminasi Gauss untuk grid 3×3 menunjukkan suhu minimum 36,5179 °C, suhu maksimum 62,2321 °C, seluruh titik berada di antara 25 °C dan 85 °C, serta hasil simetris kiri-kanan.

### 6.5 Grid 39×39 sebagai Ilustrasi

Grid 39×39 mempunyai 1.521 titik dalam dan 1.521 persamaan. Laporan menyajikannya sebagai ilustrasi distribusi suhu stabil dan grafik konvergensi, tanpa tabel angka T1 sampai T1521. Pola suhu tampak paling tinggi di dekat sisi atas dan menurun ke arah sisi bawah. Perubahan suhu tampak lebih halus karena jumlah titik dalam lebih rapat.

## 7\. Gambar Hasil Program

Menu 3 menyimpan gambar PNG ke folder `gambar\_program`.

### Gambar 1. Suhu stabil grid 3×3

!\[Suhu stabil grid 3x3](gambar\_program/gambar1\_pelat\_3x3.png)

### Gambar 2. Konvergensi grid 3×3

!\[Konvergensi grid 3x3](gambar\_program/gambar2\_konvergensi\_3x3.png)

### Gambar 3. Ilustrasi suhu grid 39×39

!\[Ilustrasi suhu grid 39x39](gambar\_program/gambar3\_peta\_suhu\_39x39.png)

### Gambar 4. Konvergensi grid 39×39

!\[Konvergensi grid 39x39](gambar\_program/gambar4\_konvergensi\_39x39.png)

## 8\. Struktur Proyek

```text
simulasi\_distribusi\_panas/
├── README.md
├── program\_chip.py
├── requirements.txt
└── gambar\_program/
    ├── gambar1\_pelat\_3x3.png
    ├── gambar2\_konvergensi\_3x3.png
    ├── gambar3\_peta\_suhu\_39x39.png
    └── gambar4\_konvergensi\_39x39.png
```

## 9\. Kesimpulan

Berdasarkan pengujian grid 3×3, ketiga metode menghasilkan suhu yang sama sampai empat angka di belakang koma. Gauss-Seidel memerlukan 27 iterasi, lebih sedikit daripada Jacobi yang memerlukan 50 iterasi untuk tebakan awal 0 °C. Eliminasi Gauss memberikan solusi langsung dengan residual sekitar (2,8\\times10^{-14}) dan digunakan sebagai solusi acuan. Tebakan awal 25 °C sedikit mengurangi jumlah iterasi tanpa mengubah hasil akhir. Waktu komputasi untuk sembilan persamaan terlalu kecil untuk dijadikan dasar utama perbandingan. Grid 39×39 digunakan sebagai ilustrasi distribusi suhu dan grafik konvergensi.

