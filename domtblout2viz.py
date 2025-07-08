import argparse
import itertools
import matplotlib.pyplot as plt

class Coord:
    from_: int
    to: int

class Hit:
    targetName: str
    eValue: float
    coord: Coord

def createArgsParser():
    parser = argparse.ArgumentParser(
                        prog='ProgramName',
                        description='What the program does',
                        epilog='Text at the bottom of help',
                        add_help=True)
    parser.add_argument('-f', '--filename', required=True, help='Input file name')
    parser.add_argument('-o', '--output', help='Output png file name. Defaults to filename-results.png')
    parser.add_argument('-l', '--limit', help='Limit the number of hits to draw. Defaults to all hits', type=int, default=None)
    return parser

def readFile(fileName, limit=None):
    data: dict[str, list[Hit]] = {}
    count = 0
    with open(fileName) as f:
        for line in f:
            if line.startswith('#'):
                continue

            parts = line.strip().split(None, 22)
            if len(parts) < 23:
                continue
            
            i_evalue = float(parts[12])
            if i_evalue > 0.05:
                continue

            target_name = parts[0]
            accession = parts[1]
            tlen = int(parts[2])
            query_name = parts[3]
            query_accession = parts[4]
            qlen = int(parts[5])
            full_seq_evalue = parts[6]
            full_seq_score = float(parts[7])
            full_seq_bias = float(parts[8])
            domain_number = int(parts[9])
            domain_total = int(parts[10])
            c_evalue = parts[11]
            domain_score = float(parts[13])
            domain_bias = float(parts[14])
            hmm_from = int(parts[15])
            hmm_to = int(parts[16])
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            env_from = int(parts[19])
            env_to = int(parts[20])
            acc = float(parts[21])
            description = parts[22]

            if query_name not in data:
                data[query_name] = []
            data[query_name].append({
                "targetName": target_name,
                "eValue": i_evalue,
                "coord": {"from_": ali_from, "to": ali_to}
            })
            count += 1
            if limit is not None and count >= limit:
                break
    return data

def filterOverlappingElements(sequenceElements:list[Hit]):
    overlappingElements: list[Hit]=[]
    noOverlappingElements: list[Hit]= []
    for i in range(len(sequenceElements)):

        currentElement= sequenceElements[i]

        currentElementCoordStart = currentElement['coord']['from_']
        currentElementCoordEnd = currentElement['coord']['to']

        overlap=False
        for j in range(i+1, len(sequenceElements)):
            currentElementToCompare= sequenceElements[j]   
            currentElementToCompareCoordStart = currentElementToCompare['coord']['from_']
            currentElementToCompareCoordEnd = currentElementToCompare['coord']['to']
            if currentElementCoordStart <= currentElementToCompareCoordEnd and currentElementToCompareCoordStart <= currentElementCoordEnd:
                overlap=True
                overlappingElements.append(currentElement)
                overlappingElements.append(currentElementToCompare)

        if not overlap and  not currentElement in overlappingElements:
            noOverlappingElements.append(currentElement)
    return {'overlappingElements':overlappingElements, 'noOverlappingElements':noOverlappingElements}

def findLowestEValueElement(overlappingElements: list[Hit]):
    lowestEValueElement:Hit=overlappingElements[0]
    for element in overlappingElements:
           if lowestEValueElement['eValue'] > element['eValue']: 
                    lowestEValueElement = element
    return lowestEValueElement

def printGraph(data, fileName):
    fig, ax = plt.subplots(figsize=(10, 6))

    line_height = 0.6  # vertical space between lines
    y_offset = 0.05    # offset for text above line

    for sequence_idx, (sequence, hits) in enumerate(data.items()):
        color_iter = itertools.cycle(plt.cm.tab10.colors)
        y = sequence_idx * line_height  # increase space between lines
        for hit_idx, hit in enumerate(hits):
            coord = hit['coord']
            color = next(color_iter)
            # Draw the line
            ax.plot([coord['from_'], coord['to']], [y, y], marker='o', linewidth=3)
            # Add the target name above the line
            mid = (coord['from_'] + coord['to']) / 2
            ax.text(mid, y + y_offset, hit['targetName'], va='bottom', ha='center', fontsize=9, color='black', fontweight='bold')
            # Add from and to coords to each line, offset vertically to avoid overlap
            ax.text(coord['from_']+5, y - 0.12, str(coord['from_']), va='top', ha='right', fontsize=8, color='gray')
            ax.text(coord['to']-5, y - 0.12, str(coord['to']), va='top', ha='left', fontsize=8, color='gray')

    yticks = [i * line_height for i in range(len(data))]
    ax.set_yticks(yticks)
    ax.set_yticklabels(list(data.keys()))
    ax.set_xlabel('Position')
    ax.set_ylabel('Sequence')
    ax.set_title('Hits Visualization')
    ax.grid(True, axis='x')
    plt.tight_layout()
    plt.savefig(fileName)
    plt.show()

if __name__ == "__main__":
    args = createArgsParser().parse_args()

    if not args.output:
        args.output=f"{args.filename}-results.png"
    print (args)

    data = readFile(args.filename, limit=args.limit)

    finalDict: dict[str, list[Hit]] = {}
    for sequence, sequenceElements in data.items():
        result= filterOverlappingElements(sequenceElements)
        overlappingElements= result['overlappingElements']
        noOverlappingElements= result['noOverlappingElements']
 
        if sequence not in finalDict:
            finalDict[sequence] = []
        finalDict[sequence].extend(noOverlappingElements)  

        lowestEValueElement:Hit
        if len(overlappingElements)>0:
            lowestEValueElement=findLowestEValueElement(overlappingElements)  
            finalDict[sequence].append(lowestEValueElement)

    printGraph(finalDict, args.output)