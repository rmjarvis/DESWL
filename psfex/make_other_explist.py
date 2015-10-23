# Make a list of run, exp that are in the full list, but not in a given partial list.

full = 'y1all'
part = 'y1spte'
other = 'y1other'
bands = ['g','r','i','z']

for b in bands:
    full_b = full + '_' + b
    part_b = part + '_' + b
    other_b = other + '_' + b

    with open(full_b,'r') as fin:
        all_runexp = set([ tuple(line.strip().split()) for line in fin ])
    with open(part_b,'r') as fin:
        some_runexp = set([ tuple(line.strip().split()) for line in fin ])

    other_runexp = all_runexp - some_runexp

    with open(other_b,'w') as fout:
        for run,exp in sorted(other_runexp):
            fout.write(run + ' ' + exp + '\n')
