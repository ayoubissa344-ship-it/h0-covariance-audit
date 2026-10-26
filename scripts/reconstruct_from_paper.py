#!/usr/bin/env python3
"""Rebuild the 88-row compilation from Appendix A (Tables X-XI) of Pantos & Perivolaropoulos,
arXiv:2601.00650v2, and write data/h0_compilation_reconstructed.csv.

This is a PROVISIONAL offline reconstruction, used when the pinned upstream script cannot be fetched.
Each row carries a `status`:
  exact        value and symmetric error as printed in Appendix A
  symmetrised  asymmetric error reduced by RMS, sigma = sqrt((s+^2 + s-^2)/2) (paper Eq. 4)
  approx       one error component was unreadable in the extracted PDF text and is estimated
               (affects only low-weight rows; see docs/reconstruction.md)
The authoritative table is the one produced by `make data` from the pinned upstream commit.
"""
import math
from pathlib import Path

rms = lambda u, d: math.sqrt((u * u + d * d) / 2.0)
# (group, label, h0, sigma, status)
T = []
def add(g, label, h0, sig, status="exact"): T.append((g, label, h0, round(sig, 6), status))
sym = lambda g, l, h, u, d: add(g, l, h, rms(u, d), "symmetrised")
apx = lambda g, l, h, u, d: add(g, l, h, rms(u, d), "approx")

# ---- Category 1 (L) : 30 rows, Table X
for l, h, s in [("Huang 2019",73.3,4.0),("de Jaeger 2020",75.8,5.1),("Kourkchi 2020",76.0,2.5),("Khetan 2021",70.50,4.13),
    ("Blakeslee 2021",73.3,2.5),("Freedman 2021",69.8,1.7),("Dhawan 2022 (a)",76.94,6.4),("Kenworthy 2022",73.1,2.5),
    ("Dhawan 2022 (b)",70.92,1.88),("Dhawan 2022 (c)",74.82,1.28),("Scolnic 2023",73.22,2.06),("de Jaeger 2023",74.1,8.0),
    ("Uddin 2023 (a)",71.76,1.32),("Uddin 2023 (b)",73.22,1.45),("Huang 2023",72.37,2.97),("Li 2024",74.7,3.1),
    ("Chavez 2024",73.1,2.3),("Boubel 2024",73.3,4.1),("Lee 2024",67.80,2.72),("Freedman 2024",72.05,3.60),
    ("Said 2024",76.05,4.90),("Zhang 2024",75.5,3.8),("Jensen 2025",73.8,2.4),("Wojtak 2025",70.59,1.15),
    ("Bhardwaj 2025",73.06,2.6),("Newman 2025",75.3,2.9),("Riess 2025",73.49,0.93),("Kudritzki 2025",76.2,6.2),
    ("Wagner 2025",68.0,8.0),("Lapuente 2025",72.61,1.69)]:
    add("L", l, h, s)
# ---- Category 2 (R) : 33 rows, Table XI
sym("R","Wu 2021",68.81,4.99,4.33); apx("R","Abbott 2021",68,8,6); sym("R","Shajib 2023",77.1,7.3,7.1); add("R","Napier 2023",74.1,8.0)
sym("R","Kelly 2023",66.6,4.1,3.3); sym("R","Dominguez 2023",61.9,2.9,2.4); apx("R","Gao 2023",65.5,6.4,5.5); sym("R","Liu 2023 (QSO cluster)",59.1,3.6,3.5)
apx("R","Fung 2023",86,55,46); sym("R","Martinez 2023",74,9,13); sym("R","Ballard 2023",85.4,29.1,33.9); add("R","Moresco 2023",66.7,5.3)
sym("R","Alfradique 2024",68.84,15.51,7.74); sym("R","Grillo 2024",66.2,3.5,3.2); sym("R","Pascale 2024",71.8,9.2,8.1)
sym("R","Hernandez 2024",73.0,13.7,7.7); sym("R","Bom 2024",70.4,13.6,11.7); sym("R","TDCOSMO XVI 2024",65,23,14); add("R","Gao 2024",68.81,5.64)
add("R","Li 2024 (lensing)",66,5); sym("R","Yang 2024",74.0,7.5,7.2); add("R","Piratova-Moreno 2025",65.1,7.4); sym("R","Barua 2025",73.5,3.0,2.9)
add("R","Zhang 2025 (FRB)",75,30); sym("R","Gao 2025 (FRB)",86.18,18.03,14.99); add("R","Loubser 2025",54.6,26.3); add("R","Beirnaert 2025",79,9)
apx("R","Bahr-Kalus 2025 (LCDM)",74.0,7.2,3.5); sym("R","Birrer 2025",71.6,3.9,3.3); add("R","Liu 2025",66.0,4.3); sym("R","Pierel 2025",66.9,11.2,8.1)
apx("R","Paic 2025",64.2,5.8,5.0); sym("R","Oliveira 2026",67.9,4.4,4.3)
# ---- Category 3 (R) : 16 rows
add("R","Kuo 2013",68,9); add("R","Gao 2015",66,6); add("R","Bulla 2022",69.6,5.5); add("R","Zhang 2022",65.9,3.0)
add("R","Du 2023",71.5,3.8); sym("R","Palmese 2023",75.46,5.34,5.39); add("R","Sneppen 2023",67.0,3.6); sym("R","Liu 2023 (lensing+SNe)",72.9,2.0,2.3)
sym("R","Gonzalez 2024",72.7,6.3,5.6); sym("R","Li 2024 (Refsdal)",66.3,3.8,3.6); apx("R","Jaiswal 2024",66.9,10.6,8.9); add("R","Vogl 2024",74.9,2.7)
sym("R","Song 2025",70.40,8.03,5.60); add("R","Colaco 2025",70.55,7.44); add("R","Du 2025",71.59,0.94); add("R","Favale 2025",68.8,3.0)
# ---- Category 4 (R) : 9 rows
add("R","Reese 2003",61,18); sym("R","Philcox 2022",64.8,2.2,2.5); add("R","Colaco 2023",67.22,6.07); add("R","Pogosian 2024 (pre-DESI)",68.05,0.94)
add("R","Pogosian 2025 (DESI DR2)",69.37,0.65); sym("R","Bahr-Kalus 2025 (CMB-NSH)",65.2,4.9,6.2); add("R","Garcia 2025",70.03,0.97)
sym("R","Zaborowski 2025",70.8,2.0,2.2); add("R","Krolewski 2025",69.0,2.5)

assert sum(r[0]=="L" for r in T)==30 and sum(r[0]=="R" for r in T)==58, (len(T))
out = Path(__file__).resolve().parents[1] / "data" / "h0_compilation_reconstructed.csv"
cnt = {k: sum(r[4]==k for r in T) for k in ("exact","symmetrised","approx")}
with open(out, "w", encoding="utf-8", newline="") as f:
    f.write("# PROVISIONAL reconstruction from Appendix A of arXiv:2601.00650v2 (see scripts/reconstruct_from_paper.py, docs/reconstruction.md).\n")
    f.write("# status: exact | symmetrised (RMS of asymmetric errors) | approx (one error component estimated). Authoritative table: `make data`.\n")
    f.write("id,group,h0,sigma,label,status\n")
    cl = cr = 0
    for g, l, h, s, st in T:
        if g == "L": cl += 1; i = f"L{cl:02d}"
        else: cr += 1; i = f"R{cr:02d}"
        f.write(f"{i},{g},{h},{s},{l.replace(',', ';')},{st}\n")
print("wrote", out, cnt)
