# 0 homozygot
# 1 heterozygot
# 2 homozygot alt
snp_marks: dict[str, dict[str, int]] = {
    "A": {"snp1": 2, "snp2": 1, "snp3": 0, "sp4": 2, "snp5": 1},
    "B": {"snp1": 2, "snp2": 1, "snp3": 2, "sp4": 2, "snp5": 1},
    "C": {"snp1": 0, "snp2": 2, "snp3": 1, "sp4": 0, "snp5": 1}
}
#find pairs of individuals that's most like to each other based on the snp_marks

individuals = list(snp_marks.keys())

# marks shared:
#A & B: snp1, snp2, sp4, snp5
#A & C: snp2, sp4, snp5
#B & C: snp2, sp4, snp5
def find_pairs(individuals: list[str])-> set[tuple[str, str]]:
    pairs: set[tuple[str, str]] = set()

    for i in range(len(individuals)):
        for j in range(i + 1, len(individuals)):
            pairs.add((individuals[i], individuals[j]))

    print(pairs)
    return pairs
def shared_snp_count(ind1: str, ind2: str) -> int:
    shared_count = 0
    for snp in snp_marks[ind1]:
        #for i in range(len(snp_marks[ind1])):
        if snp_marks[ind1][snp] == snp_marks[ind2][snp]:
            shared_count += 1
    return shared_count
def most_similar_pairs(pairs: set[tuple[str, str]]) -> list[tuple[str, float]]:
    similarity_scores: list[tuple[str, float]] = []
    for ind1, ind2 in pairs:
        shared_count = shared_snp_count(ind1, ind2)
        print(f"Shared SNPs between {ind1} and {ind2}: {shared_count}")
        total_snps = len(snp_marks[ind1])
        similarity = shared_count / total_snps
        similarity_scores.append((f"{ind1} & {ind2}", similarity))
    return sorted(similarity_scores, key=lambda x: x[1], reverse=True)

def main():
    pairs = find_pairs(individuals)
    most_similar = most_similar_pairs(pairs)
    result = max(most_similar, key=lambda x: x[1])
    print(f"Most similar pair: {result[0]}, Similarity Score: {result[1]:.2f}")

if __name__ == "__main__":
    main()