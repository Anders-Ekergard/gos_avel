"""
OCS-optimering (steg N i pipelineskissen) - forenklad Optimal Contribution
Selection. Kombinerar avelsvarde (steg M) och genomisk likhet/GRM (steg L)
for att foresla avelspar: maximera genetisk vinst men undvik att para de
mest slaktskapsnara individerna.

OBS: leksaksdata (n=20). En riktig OCS (Meuwissen 1997) loser ett
kontinuerligt optimeringsproblem (kvadratisk programmering) over hur mycket
varje individ ska bidra till nasta generation, under bivillkoret att
medelslaktskapet halls under ett mal. Har anvands en girig approximation
som visar samma PRINCIP i mindre skala: valj de basta individerna (lagst
avelsvarde), para dem tva och tva sa att inget par overskrider en
slaktskapstroskel.

Kor: python ocs.py
"""

from itertools import combinations

from marker_fenotyp_koppling import GENOTYP_CSV, FENOTYP_CSV, las_genotyper, las_fenotyper
from avelsvarde import markor_effekter, berakna_avelsvarden
from gos_genomik_resurser import (
    INAVELSKOEFFICIENT_F_AVKOMMA as PUBLICERAD_F,
    EFFEKTIV_POPULATIONSSTORLEK_NE as PUBLICERAD_NE,
    SRA_PROJECT_HALVSYSKON,
)

# Framtida valideringsmal for den har modulen (inte implementerat an): en
# publicerad studie pa riktig gos-halvsyskondata (SRA-projekt, se
# SRA_PROJECT_HALVSYSKON ovan) fann F=0.33, Ne=12 (se PUBLICERAD_F/
# PUBLICERAD_NE ovan) - allvarligt lagt, langt under rekommenderat
# minimum ~50. Om/nar genotyper/fenotyper fran den datan gors
# tillgangliga har (kraver egen alignment/variant-calling-pipeline pa
# rå-FASTQ, se gos_genomik_resurser.py), ar ett konkret sanity-check att
# kora denna moduls f_prognos/genomisk_likhet-logik mot den riktiga
# datan och se om den ateraper samma storleksordning - inte bara mot
# leksaksdatan (n=20) som anvands nu.

# Hur manga av de basta individerna (efter avelsvarde) som blir avelskandidater
ANTAL_KANDIDATER = 8

# Andel delade alleldoser (identity-by-state) over vilken ett par undviks,
# oavsett hur bra deras avelsvarden ar
SLAKTSKAP_TROSKEL = 0.7

# Hur manga generationer F-prognosen projicerar framat
GENERATIONER = 5


def genomisk_likhet(genotyper: dict[str, dict[str, str]], ind1: str, ind2: str) -> float:
    """Andel markorer dar tva individer har identisk genotyp - samma
    princip som grm.py, men over alla markorer i CSV-datan istallet for
    det hardkodade A/B/C-exemplet."""
    markorer = genotyper[ind1].keys()
    delade = sum(1 for m in markorer if genotyper[ind1][m] == genotyper[ind2][m])
    return delade / len(markorer)


def foresla_par(
    kandidater: list[str],
    genotyper: dict[str, dict[str, str]],
    avelsvarden: dict[str, float],
    slaktskap_troskel: float = SLAKTSKAP_TROSKEL,
) -> list[tuple[str, str, float]]:
    """Girig parning: ta individerna i avelsvarde-ordning, para varje
    obunden individ med basta obundna individ som INTE ligger over
    slaktskapstroskeln. Finns ingen sadan, valj minst slaktskapsnara."""
    aterstar = sorted(kandidater, key=lambda d: avelsvarden[d])
    par: list[tuple[str, str, float]] = []

    while len(aterstar) >= 2:
        forsta = aterstar.pop(0)
        basta_partner = None
        for kandidat in aterstar:
            if genomisk_likhet(genotyper, forsta, kandidat) <= slaktskap_troskel:
                basta_partner = kandidat
                break
        if basta_partner is None:
            basta_partner = min(aterstar, key=lambda k: genomisk_likhet(genotyper, forsta, k))
        aterstar.remove(basta_partner)
        likhet = genomisk_likhet(genotyper, forsta, basta_partner)
        par.append((forsta, basta_partner, likhet))

    return par


def naiv_parning(kandidater: list[str], avelsvarden: dict[str, float]) -> list[tuple[str, str]]:
    """Jamforelse-baseline utan slaktskapshansyn: sortera efter avelsvarde
    och para sekventiellt (basta med nastbasta, osv). Visar vad OCS-steget
    faktiskt tillfor jamfort med att bara jaga avelsvarde."""
    ordnade = sorted(kandidater, key=lambda d: avelsvarden[d])
    return [(ordnade[i], ordnade[i + 1]) for i in range(0, len(ordnade) - 1, 2)]


def f_prognos(medel_f_idag: float, generationer: int, ne: int) -> list[float]:
    """Standardapproximation fran kvantitativ genetik: forandringstakt i
    inavel per generation ~ 1/(2*Ne). F vaxer da enligt
    F_t = F_(t-1) + (1 - F_(t-1)) * delta_F.

    OBS: Ne har skattas grovt som antalet avelsdjur som faktiskt anvands
    (ingen konnsfordelning i leksaksdatan) - en verklig Ne blir ofta lagre
    an sa (skev konnsfordelning/familjestorlek sanker den). Prognosen
    antar ocksa att samma strategi och samma Ne aterupprepas varje
    generation - en trend, inte en simulering av faktiska framtida par."""
    delta_f = 1 / (2 * ne)
    prognos = []
    f = medel_f_idag
    for _ in range(generationer):
        f = f + (1 - f) * delta_f
        prognos.append(f)
    return prognos


def kor_avelspipeline(
    genotyper: dict[str, dict[str, str]],
    fenotyper: dict[str, float],
    antal_kandidater: int = ANTAL_KANDIDATER,
    slaktskap_troskel: float = SLAKTSKAP_TROSKEL,
    generationer: int = GENERATIONER,
) -> dict:
    """Kor hela Del 2-kedjan (avelsvarde -> GRM/genomisk likhet ->
    OCS-parforslag -> F-prognos) och returnerar strukturerad data, utan
    att skriva ut nagot.

    Ren datainsamlingsfunktion - delad av CLI:t (main) och webb-API:t,
    samma monster som candidate_gene_dossier.gather_dossier, sa de tva
    ytorna aldrig kan divergera i vad de faktiskt rapporterar."""
    effekter = markor_effekter(genotyper, fenotyper)
    varden = berakna_avelsvarden(genotyper, effekter)

    kandidater = sorted(varden, key=lambda d: varden[d])[:antal_kandidater]

    likhetsmatris = {
        f"{a}|{b}": genomisk_likhet(genotyper, a, b)
        for a, b in combinations(kandidater, 2)
    }

    par = foresla_par(kandidater, genotyper, varden, slaktskap_troskel)
    ocs_medel_f = sum(likhet for _, _, likhet in par) / len(par)

    naiva_par = naiv_parning(kandidater, varden)
    naiv_medel_f = sum(genomisk_likhet(genotyper, a, b) for a, b in naiva_par) / len(naiva_par)

    ne = len(par) * 2
    prognos = f_prognos(ocs_medel_f, generationer, ne)

    return {
        "markor_effekter": {
            m: {"effekt": e["effekt"], "additiv_avvikelse": e["additiv_avvikelse"]}
            for m, e in effekter.items()
        },
        "avelsvarden": [
            {"djur": d, "avelsvarde": varden[d]} for d in sorted(varden, key=lambda d: varden[d])
        ],
        "kandidater": kandidater,
        "likhetsmatris": likhetsmatris,
        "slaktskap_troskel": slaktskap_troskel,
        "foreslagna_par": [
            {
                "a": a,
                "b": b,
                "f_proxy": likhet,
                "over_troskel": likhet > slaktskap_troskel,
                "avelsvarde_a": varden[a],
                "avelsvarde_b": varden[b],
            }
            for a, b, likhet in par
        ],
        "ocs_medel_f": ocs_medel_f,
        "naiv_medel_f": naiv_medel_f,
        "ne": ne,
        "generationer": generationer,
        "f_prognos": prognos,
        "publicerad_referens": {
            "f": PUBLICERAD_F,
            "ne": PUBLICERAD_NE,
            "sra_projekt": SRA_PROJECT_HALVSYSKON,
        },
    }


def main():
    genotyper = las_genotyper(GENOTYP_CSV)
    fenotyper = las_fenotyper(FENOTYP_CSV)

    resultat = kor_avelspipeline(genotyper, fenotyper)

    print(f"--- {ANTAL_KANDIDATER} avelskandidater (lägst avelsvärde = bäst) ---")
    for djur in resultat["kandidater"]:
        varde = next(a["avelsvarde"] for a in resultat["avelsvarden"] if a["djur"] == djur)
        print(f"  {djur}: avelsvärde={varde:+.3f}")
    print()

    print(f"--- Genomisk likhet mellan kandidaterna (tröskel: {SLAKTSKAP_TROSKEL}) ---")
    for par_namn, likhet in resultat["likhetsmatris"].items():
        a, b = par_namn.split("|")
        flagga = "  <- över tröskeln" if likhet > SLAKTSKAP_TROSKEL else ""
        print(f"  {a} & {b}: {likhet:.2f}{flagga}")
    print()

    print("--- Föreslagna avelspar (OCS) ---")
    for p in resultat["foreslagna_par"]:
        varning = ""
        if p["over_troskel"]:
            varning = "  <- ingen partner under tröskeln fanns bland återstående, minst släktskapsnära valdes"
        print(
            f"  {p['a']} x {p['b']}  (F-proxy avkomma: {p['f_proxy']:.2f}, "
            f"avelsvärden: {p['avelsvarde_a']:+.3f} / {p['avelsvarde_b']:+.3f}){varning}"
        )
    print()

    print("--- F-jämförelse: OCS vs naiv parning (bara efter avelsvärde) ---")
    print(f"  OCS medel-F denna generation:  {resultat['ocs_medel_f']:.3f}")
    print(f"  Naiv medel-F denna generation: {resultat['naiv_medel_f']:.3f}")
    if resultat["naiv_medel_f"] > 0:
        minskning = (resultat["naiv_medel_f"] - resultat["ocs_medel_f"]) / resultat["naiv_medel_f"]
        print(f"  -> OCS sänker snittet med {minskning:.0%}")
    print()

    print(f"--- F-prognos, OCS-strategin upprepad (Ne~{resultat['ne']}) ---")
    for generation, f in enumerate(resultat["f_prognos"], start=1):
        print(f"  Generation +{generation}: F~{f:.3f}")
    print()
    print("OBS: grov trend (konstant Ne, ingen faktisk simulering av framtida")
    print("genotyper) - se docstring för f_prognos för antaganden.")
    print()
    ref = resultat["publicerad_referens"]
    print(f"Jämförelsepunkt (riktig data, ej denna leksaksdata): en publicerad studie på")
    print(f"gös-halvsyskon (SRA {ref['sra_projekt']}) fann F={ref['f']}, Ne={ref['ne']}.")


if __name__ == "__main__":
    main()
