"""
Avelsvarde per individ (steg M i pipelineskissen) - kombinerar markorernas
additiva effekter (steg K, marker_fenotyp_koppling.py) till ett sammanvagt
varde per djur.

OBS: leksaksdata (n=20) - visar PRINCIPEN for hur ett avelsvarde racknas
fram fran markoreffekter, ingen slutsats om gos.

Kor: python avelsvarde.py
"""

from marker_fenotyp_koppling import GENOTYP_CSV, FENOTYP_CSV, las_genotyper, las_fenotyper, medel

# Alleldos: antal "B"-alleler vid markoren (0 = AA, 1 = AB, 2 = BB)
ALLELDOS = {"AA": 0, "AB": 1, "BB": 2}


def markor_effekter(genotyper: dict[str, dict[str, str]], fenotyper: dict[str, float]) -> dict[str, dict]:
    """For varje markor: additiv effekt per "B"-allel (halva skillnaden
    mellan homozygotgrupperna), plus hur val heterozygoten foljer en
    additiv modell (avstand fran forvantad mittpunkt (AA+BB)/2)."""
    djur_ids = [d for d in genotyper if d in fenotyper]
    markor_namn = next(iter(genotyper.values())).keys()

    effekter = {}
    for markor in markor_namn:
        grupper: dict[str, list[float]] = {"AA": [], "AB": [], "BB": []}
        for djur in djur_ids:
            grupper[genotyper[djur][markor]].append(fenotyper[djur])

        if not grupper["AA"] or not grupper["BB"]:
            effekter[markor] = {"effekt": 0.0, "additiv_avvikelse": None}
            continue

        medel_aa = medel(grupper["AA"])
        medel_bb = medel(grupper["BB"])
        effekt = (medel_bb - medel_aa) / 2  # forandring i FCR per "B"-allel

        additiv_avvikelse = None
        if grupper["AB"]:
            forvantad_ab = (medel_aa + medel_bb) / 2
            additiv_avvikelse = medel(grupper["AB"]) - forvantad_ab

        effekter[markor] = {"effekt": effekt, "additiv_avvikelse": additiv_avvikelse}
    return effekter


def berakna_avelsvarden(
    genotyper: dict[str, dict[str, str]], effekter: dict[str, dict]
) -> dict[str, float]:
    """Per individ: summan av alleldos * markoreffekt over alla markorer.
    Lagre varde = lagre forvantad FCR = battre, eftersom foderkvot ska
    minimeras."""
    varden = {}
    for djur, dgeno in genotyper.items():
        varde = 0.0
        for markor, genotyp in dgeno.items():
            dos = ALLELDOS[genotyp]
            varde += dos * effekter[markor]["effekt"]
        varden[djur] = varde
    return varden


def main():
    genotyper = las_genotyper(GENOTYP_CSV)
    fenotyper = las_fenotyper(FENOTYP_CSV)

    effekter = markor_effekter(genotyper, fenotyper)

    print("--- Markoreffekter (additiv modell: AA vs BB) ---")
    for markor, e in effekter.items():
        varning = ""
        if e["additiv_avvikelse"] is not None and abs(e["additiv_avvikelse"]) > 0.15:
            varning = "  <- heterozygot följer INTE additiv modell väl, använd med försiktighet"
        print(f"  {markor}: effekt/B-allel={e['effekt']:+.3f}{varning}")
    print()

    varden = berakna_avelsvarden(genotyper, effekter)

    print("--- Avelsvärde per djur (lägre = bättre, förväntat FCR-bidrag) ---")
    for djur, varde in sorted(varden.items(), key=lambda x: x[1]):
        print(f"  {djur}: avelsvärde={varde:+.3f}  (faktisk FCR: {fenotyper[djur]:.2f})")

    print()
    print("OBS: alla markörer räknas in, även de vars heterozygot inte följer")
    print("en additiv modell väl (flaggade ovan) - deras bidrag är därför brus,")
    print("inte signal. I ett riktigt urval bör sådana markörer viktas ner eller")
    print("uteslutas tills mer data bekräftar dem.")


if __name__ == "__main__":
    main()
