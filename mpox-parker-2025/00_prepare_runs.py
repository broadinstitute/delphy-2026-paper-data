#!/usr/bin/env python3

from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET
import zipfile

# Prepare directory structure
# ===========================
raw_path = Path('raw')
raw_path.mkdir(parents=True, exist_ok=True)
delphy_inputs_path = Path('delphy_inputs')
delphy_inputs_path.mkdir(parents=True, exist_ok=True)

# Clone O'Toole et al data repo
# =============================
print("\nCloning Parker et al 2025 data repo")
data_repo_path = raw_path / 'Mpox_West_Africa'
do_clone = True
if data_repo_path.exists():
    do_clone = False
    print(f'[SKIPPING] {data_repo_path} already exists')

if do_clone:
    subprocess.run([
        'git',
        'clone',
        'https://github.com/andersen-lab/Mpox_West_Africa.git'],
                   cwd=raw_path.as_posix())

orig_beast_zip_file_path = data_repo_path / 'BEAST' / 'Mpox_2epoch_combinedDTA.xml.zip'
if not orig_beast_zip_file_path.exists():
    raise ValueError(f"Zipped BEAST XML file not found in Parker et al 2025 data repo: {orig_beast_zip_file_path}")

# Extract sequences, states and regions from BEAST file
# =====================================================
print("\nReading sequence data from BEAST XML file...")
with zipfile.ZipFile(orig_beast_zip_file_path.as_posix(), 'r') as archive:
    with archive.open('Mpox_2epoch_combinedDTA.xml', 'r') as f:
        beastXml = ET.parse(f)
alignments = beastXml.findall("alignment")
apobec3_alignment = [a for a in alignments if a.get('id') == 'apobec3_alignment'][0]
non_apobec3_alignment = [a for a in alignments if a.get('id') == 'non_apobec3_alignment'][0]

included_sequences = []

# Check that sequences are listed in the same order in both alignments
for apo, nonapo in zip(apobec3_alignment.findall('sequence'),
                       non_apobec3_alignment.findall('sequence')):
    apo_id = apo.find("taxon").get("idref")
    non_apo_id = nonapo.find("taxon").get("idref")
    
    if apo_id != non_apo_id:
        raise ValueError(f'Sequences in BEAST files alignments listed in different order (first discrepancy: {apo_id} vs {non_apo_id})')


# First extract state & region from BEAST XML
id_2_state = {}
id_2_region = {}

taxas = beastXml.findall("taxa")
main_taxa_xml = [t for t in taxas if t.get('id') == 'taxa'][0]  # First taxa element with id "taxa"
for taxon_xml in main_taxa_xml:
    theId = taxon_xml.get('id')
    attrs = taxon_xml.findall('attr')
    id_2_state[theId] = [a for a in attrs if a.get('name') == 'State'][0].text.strip()
    id_2_region[theId] = [a for a in attrs if a.get('name') == 'Region'][0].text.strip()

if False:   # Debugging
    for state in sorted(set(id_2_state.values())):
        print(f'State: {state}')
        for theId, theState in id_2_state.items():
            if state == theState:
                print(f'- {theId}')
    
    for region in sorted(set(id_2_region.values())):
        print(f'Region: {region}')
        for theId, theRegion in id_2_region.items():
            if region == theRegion:
                print(f'- {theId}')

# Exclude sequences not in scope (GISAID, pre-spillover, non-CladeIIb)
def validate_date(date_str):
    if len(date_str) == 4:         # Samples without exact month & day...
        return date_str  # No problem now that we have tip-date sampling
    #    return date_str + '-07-01' # ...are dated to mid-year
    if len(date_str) == 4+1+2:     # Samples without exact day...
        return date_str  # No problem now that we have tip-date sampling
    #    return date_str + '-15'    # ...are dated to mid-month
    if len(date_str) != 4+1+2+1+2:
        raise ValueError(f'Invalid date string: {date_str}')
    return date_str

for d in apobec3_alignment.findall('sequence'):
    theId = d.find("taxon").get('idref')
    theDate = validate_date(theId.split('|')[-1])  # All dates in the BEAST XML pass, so we use them as-is later

    # Skip GISAID sequences
    if theId.startswith('EPI'):
        print(f'Skipping GISAID sequence {theId}')
        continue

    # Skip pre-spillover samples, including any sequences before 2017
    # (e.g., outgroups KJ642615 from 1978 and KJ642617 from 1971)
    yy = int(theDate[0:4])
    if yy < 2017 or any(x in theId for x in ['VSP189', 'VSP191', 'VSP199', 'TRM288']):
        print(f'Skipping pre-spillover sequence {theId}')
        continue

    # Skip non-Clade-IIb samples
    if any(x in theId for x in [
            'TRM054',
            'TRM055',
            'TRM057',
            'TRM061',
            'TRM063',
            'TRM064',
            'TRM072',
            'TRM073',
            'TRM074',
    ]):
        print(f'Skipping non-Clade-IIb sequence {theId}')
        continue

    print(f'Found sequence for {theId}, dated to {theDate}')
    
    included_sequences.append(theId)
    
# Prepare fastas with all the sequences
# =====================================
print("\nReconstituting multiple sequence alignment from 2 partitions in BEAST XML file...")

def shorten_id(long_id):
    if long_id.split('|')[0] == 'unpub':   # e.g., unpub|TRM076|Nigeria|Rivers|2022-09-27
        return long_id.split('|')[1]
    elif long_id.startswith('PP'):         # e.g., PP852976|MPXV|TRM081|Nigeria|Bayelsa|2022-10-21
        return long_id.split('|')[2]
    else:
        return long_id.split('|')[0]  # e.g., OP535319|MPXV|Nigeria|Cross-River|2017-12

input_fasta_path = delphy_inputs_path / 'mpox-parker-2025.fasta'
do_reassembly = True
if input_fasta_path.exists():
    do_reassembly = False
    print(f'[SKIPPING] {input_fasta_path} already exists')

if do_reassembly:
    with open(input_fasta_path, 'w') as f:
        
        for apo, nonApo in zip(apobec3_alignment.findall('sequence'),
                               non_apobec3_alignment.findall('sequence')):
            apoId = apo.find("taxon").get("idref")
            nonApoId = nonApo.find("taxon").get("idref")
            
            if apoId != nonApoId:
                raise ValueError('APO and non-APO sequences with different ids?')
            
            lastBarIndex = apoId.rindex('|')
            theId, theDate = apoId[:lastBarIndex], apoId[lastBarIndex+1:]
    
            if apoId not in included_sequences:
                continue

            apoSeq = ''.join(t for t in apo.itertext()).strip()
            nonApoSeq = ''.join(t for t in nonApo.itertext()).strip()
            states = []
            for apoState, nonApoState in zip(apoSeq, nonApoSeq):
                if apoState == 'N':
                    states.append(nonApoState)
                else:
                    states.append(apoState)
            
            seq = ''.join(states)
                
            f.write('>')
            fastaId = f'{shorten_id(apoId)}|{theId}|{theDate}'  # Extract simple unique id for web UI and metadata
            f.write(fastaId)
            f.write('\n')
            f.write(seq)
            f.write('\n')
                

# Prepare metadata file
# =====================
print("\nPreparing metadata...")
input_metadata_path = delphy_inputs_path / 'mpox-parker-2025_metadata.csv'
do_metadata = True
if input_metadata_path.exists():
    do_metadata = False
    print(f'[SKIPPING] {input_metadata_path} already exist')

if False:  # Debugging
    print(f'States: {sorted(set(id_2_state.values()))}')
    print(f'Regions: {sorted(set(id_2_region.values()))}')

state_2_clean_state = {
    'Akwa-Ibom':      'Akwa-Ibom',
    'Cross-River':    'Cross-River',
    'abia':           'Abia',
    'anambra':        'Anambra',
    'bayelsa':        'Bayelsa',
    'benue':          'Benue',
    'borno':          'Borno',
    'delta':          'Delta',
    'ebonyi':         'Ebonyi',
    'edo':            'Edo',
    'enugu':          'Enugu',
    'fct':            'FCT',
    'gombe':          'Gombe',
    'imo':            'Imo',
    'kaduna':         'Kaduna',
    'kano':           'Kano',
    'kogi':           'Kogi',
    'kwara':          'Kwara',
    'lagos':          'Lagos',
    'nasarawa':       'Nasarawa',
    'niger':          'Niger',
    'ogun':           'Ogun',
    'ondo':           'Ondo',
    'oyo':            'Oyo',
    'plateau':        'Plateau',
    'rivers':         'Rivers',
    
    'nigeria':        'Nigeria',
    
    'cameroon':       'Rest',
    'usa':            'Rest',
}

region_2_clean_region = {
    'NC':             'NC',
    'NE':             'NE',
    'NW':             'NW',
    'S':              'SS',
    'SE':             'SE',
    'SW':             'SW',
    
    'nigeria':        'Nigeria',
    
    'Rest':           'Rest',
    'cameroon':       'Rest',
}

if do_metadata:
    with open(input_metadata_path, 'w') as f:
        f.write('id,Region,State,rawFastaId\n')
        for theId in included_sequences:
            seqId = shorten_id(theId)
            region = region_2_clean_region[id_2_region[theId]]
            state = state_2_clean_state[id_2_state[theId]]
            
            f.write(f'{seqId},{region},{state},{theId}\n')

    print(f'Metadata written to {input_metadata_path}')

# Prepare "equivalent" BEAST run for just the hMPXV-1 sequences (no DTA, no pre-spillover samples, no transition)
# ===============================================================================================================
# This is a bit nasty, but at least you see all the changes clearly.
#
# WARNING: beastXml is modified in place, so won't reflect the original contents of the BEAST XML
# from the mpox paper repo after this block.   That's why this is the last step.

beast_run_path = Path('beastX_run')
beast_run_path.mkdir(parents=True, exist_ok=True)

print("\nPreparing \"equivalent\" hMPXV-1 BEAST X 10.5.0 XML file...")
beast_path = beast_run_path / 'mpox-parker-2025-beastX.xml'
do_beast = True
if beast_path.exists():
    do_beast = False
    print(f'[SKIPPING] {beast_path} already exists')

if do_beast:
    root = beastXml.getroot()

    # Remove excluded taxa
    for taxa in root.findall('taxa'):
        if taxa.get('id') != 'taxa':
            continue
        to_remove = [t for t in taxa.findall('taxon') if t.get('id') not in included_sequences]
        for t in to_remove:
            #print(f'Removing {t.get("id")}')
            taxa.remove(t)
            
    # Remove excluded taxa from alignments
    for alignment in root.findall('alignment'):
        to_remove = []
        for sequence in alignment.findall('sequence'):
            seq_id = sequence.find('taxon').get('idref')
            if seq_id not in included_sequences:
                #print(f'Removing {seq_id} from {alignment.get("id")}')
                to_remove.append(sequence)
        for s in to_remove:
            alignment.remove(s)
        
    # Remove excluded taxa from treeModel's uncertain leafHeights
    treeModel = root.find('treeModel')
    to_remove = []
    for elem in treeModel:
        if elem.tag == 'leafHeight' and elem.get('taxon') not in included_sequences:
            to_remove.append(elem)
    
    for e in to_remove:
        treeModel.remove(e)

    # Remove DTA attributes from taxa
    for taxa in root.findall('taxa'):
        if taxa.get('id') != 'taxa':
            continue
        for taxon in taxa.findall('taxon'):
            subelems = list(taxon.iter())
            for elem in subelems:
                if elem.tag == 'attr':
                    taxon.remove(elem)
                    
    # Remove top-level "ingroup" and "outgroup" taxa and everything associated with them
    to_remove = [t for t in beastXml.findall('taxa') if t.get('id') in ["ingroup", "outgroup"]]
    for t in to_remove:
        root.remove(t)
    
    to_remove = []
    for tmrcaStatistic in root.findall('tmrcaStatistic'):
        theId = tmrcaStatistic.get('id')
        if 'ingroup' in theId or 'outgroup' in theId:
            to_remove.append(tmrcaStatistic)
    
    for coal in root.findall('coalescentLikelihood'):
        if coal.get('id') == 'coalescent.outgroup':
            to_remove.append(coal)
    
        if coal.get('id') == 'coalescent.ingroup':
            coal.set('id', 'coalescent')
            coal.find('include').find('taxa').set('idref', 'taxa')
    
    for m in root.findall('localClockModel'):
        if m.get('id') == 'apobec3.branchRates':
            m.clear()
            m.tag = 'strictClockBranchRates'
            m.set('id', 'apobec3.branchRates')
            r = ET.Element('rate')
            r.append(ET.Element('parameter', {'idref': 'apobec3.clock.rate'}))
            m.append(r)
    
    for e in root.findall('productParameter'):
        if e.get('id') == 'apobec3.stem.time':
            to_remove.append(e)
    
    for e in root.findall('sumParameter'):
        if e.get('id') == 'apobec3.transition.time':
            to_remove.append(e)
    
    for e in root.findall('ageStatistic'):
        if e.get('id') == 'age(apobec3.transition)':
            to_remove.append(e)
    
    for e in root.findall('rateStatistic'):
        if e.get('id') == 'apobec3.meanRate':
            e.find('localClockModel').tag = 'strictClockBranchRates'
    for e in root.findall('ancestralTreeLikelihood'):
        e.tag = 'treeDataLikelihood'
        if e.get('id') == 'apobec3.treeLikelihood':
            e.find('localClockModel').tag = 'strictClockBranchRates'
    
    for e in to_remove:
        root.remove(e)

    # Remove DTA-related elements
    to_remove = []
    for e in root.findall('generalDataType'):
        if e.get('id') in ['Region.dataType', 'State.dataType']:
            to_remove.append(e)
    for e in root.findall('attributePatterns'):
        if e.get('id') in ['Region.pattern', 'State.pattern']:
            to_remove.append(e)
    for e in root.findall('generalSubstitutionModel'):
        if e.get('id') in ['Region.model', 'State.model']:
            to_remove.append(e)
    for e in root.findall('sumStatistic'):
        if e.get('id') in ['Region.nonZeroRates', 'State.nonZeroRates']:
            to_remove.append(e)
    for e in root.findall('productStatistic'):
        if e.get('id') in ['Region.actualRates', 'State.actualRates']:
            to_remove.append(e)
    for e in root.findall('siteModel'):
        if e.get('id') in ['Region.siteModel', 'State.siteModel']:
            to_remove.append(e)
    for e in root.findall('markovJumpsTreeLikelihood'):
        if e.get('id') in ['Region.treeLikelihood', 'State.treeLikelihood']:
            to_remove.append(e)
    
    for e in to_remove:
        root.remove(e)
    
    
    
    mcmc = root.find('mcmc')
    joint = mcmc.find('joint')
    
    prior = joint.find('prior')
    to_remove = []
    for elem in prior:
        # DTA
        if (elem.tag == 'ctmcScalePrior' and 
            elem.find('ctmcScale').find('parameter').get('idref') in ['Region.clock.rate', 'State.clock.rate']):
            to_remove.append(elem)
        if elem.tag == 'poissonPrior' and elem.find('statistic').get('idref') in ['Region.nonZeroRates', 'State.nonZeroRates']:
            to_remove.append(elem)
        if elem.tag == 'uniformPrior' and elem.find('parameter').get('idref') in ['Region.frequencies', 'State.frequencies']:
            to_remove.append(elem)
        if elem.tag == 'cachedPrior' and elem.find('parameter').get('idref') in ['Region.rates', 'State.rates']:
            to_remove.append(elem)
        if elem.tag == 'uniformPrior' and elem.find('parameter').get('idref') in ['Region.root.frequencies', 'State.root.frequencies']:
            to_remove.append(elem)
        if elem.tag == 'strictClockBranchRates' and elem.get('idref') in ['Region.branchRates', 'State.branchRates']:
            to_remove.append(elem)
        if elem.tag == 'generalSubstitutionModel' and elem.get('idref') in ['Region.model', 'State.model']:
            to_remove.append(elem)

        # Spillover
        if elem.tag == 'uniformPrior' and elem.find('parameter').get('idref') == 'apobec3.stem.proportion':
            to_remove.append(elem)
        if elem.tag == 'oneOnXPrior' and elem.find('parameter').get('idref') == 'constant.popSize':
            to_remove.append(elem)
        if elem.tag == 'coalescentLikelihood' and elem.get('idref') == 'coalescent.outgroup':
            to_remove.append(elem)
        if elem.tag == 'coalescentLikelihood' and elem.get('idref') == 'coalescent.ingroup':
            elem.set('idref', 'coalescent')
        if elem.tag == 'localClockModel' and elem.get('idref') == 'apobec3.branchRates':
            elem.tag = 'strictClockBranchRates'
        
    for e in to_remove:
        prior.remove(e)
    
    
    likelihood = joint.find('likelihood')
    to_remove = []
    for elem in likelihood:
        # DTA
        if elem.tag == 'markovJumpsTreeLikelihood' and elem.get('idref') in ['Region.treeLikelihood', 'State.treeLikelihood']:
            to_remove.append(elem)
        if elem.tag == 'ancestralTreeLikelihood' and elem.get('idref') == 'apobec3.treeLikelihood':
            elem.tag = 'treeDataLikelihood'
        
    for e in to_remove:
        likelihood.remove(e)

    to_remove = []
    for elem in root:
        if elem.tag == 'strictClockBranchRates' and elem.get('id') in ['Region.branchRates', 'State.branchRates']:
            to_remove.append(elem)
        if elem.tag == 'rateStatistic' and elem.get('id') in ['Region.meanRate', 'State.meanRate']:
            to_remove.append(elem)
    
    for t in to_remove:
        root.remove(t)
    
    operators = root.find('operators')
    to_remove = []
    for op in operators:
        # DTA
        if op.tag == 'scaleOperator' and op.find('parameter').get('idref') in ['Region.clock.rate', 'State.clock.rate']:
            to_remove.append(op)
        if op.tag == 'upDownOperator' and op.find('down').find('parameter').get('idref') in ['Region.clock.rate', 'State.clock.rate']:
            to_remove.append(op)
        if op.tag == 'scaleOperator' and op.find('parameter').get('idref') in ['Region.rates', 'State.rates']:
            to_remove.append(op)
        if op.tag == 'bitFlipOperator' and op.find('parameter').get('idref') in ['Region.indicators', 'State.indicators']:
            to_remove.append(op)
        if op.tag == 'deltaExchange' and op.find('parameter').get('idref') in ['Region.root.frequencies', 'State.root.frequencies']:
            to_remove.append(op)
    
        # Spillover
        if op.tag == 'scaleOperator' and op.find('parameter').get('idref') == 'constant.popSize':
            to_remove.append(op)
        if op.tag == 'randomWalkOperator' and op.find('parameter').get('idref') == 'apobec3.stem.proportion':
            to_remove.append(op)
        if op.tag == 'uniformOperator' and op.find('parameter').get('idref') not in [f'age({s})' for s in included_sequences]:
            to_remove.append(op)

    for t in to_remove:
        operators.remove(t)

    to_remove = []
    for log in mcmc.findall('log'):
        if log.get('id') in ['Mpox_2poch_combined.RegionrateMatrixLog', 'Mpox_2poch_combined.StaterateMatrixLog']:
            to_remove.append(log)
    for logTree in mcmc.findall('logTree'):
        if logTree.get('fileName') in ['Mpox_2poch_combined.Region.history.trees', 'Mpox_2poch_combined.State.history.trees']:
            to_remove.append(logTree)
        if logTree.get('id') == 'treeFileLog':
            for trait in logTree.findall('trait'):
                if trait.get('tag') == 'apobec3.rate':
                    for m in trait.findall('localClockModel'):
                        m.tag = 'strictClockBranchRates'
    for l in to_remove:
        mcmc.remove(l)
    
    for log in mcmc.findall('log'):
        to_remove = []
        for elem in log:
            # DTA
            if elem.tag == 'column' and elem.get('label') in ['Region.clock.rate', 'State.clock.rate']:
                to_remove.append(elem)
            if elem.tag == 'column' and elem.get('label') in ['Region.nonZeroRates', 'State.nonZeroRates']:
                to_remove.append(elem)
            if elem.tag == 'rateStatistic' and elem.get('idref') in ['Region.meanRate', 'State.meanRate']:
                to_remove.append(elem)
            if elem.tag == 'parameter' and elem.get('idref') in ['Region.rates', 'State.rates']:
                to_remove.append(elem)
            if elem.tag == 'parameter' and elem.get('idref') in ['Region.indicators', 'State.indicators']:
                to_remove.append(elem)
            if elem.tag == 'parameter' and elem.get('idref') in ['Region.nonZeroRates', 'State.nonZeroRates']:
                to_remove.append(elem)
            if elem.tag == 'parameter' and elem.get('idref') in ['Region.clock.rate', 'State.clock.rate']:
                to_remove.append(elem)
            if elem.tag == 'strictClockBranchRates' and elem.get('idref') in ['Region.branchRates', 'State.branchRates']:
                to_remove.append(elem)
            if elem.tag == 'sumStatistic' and elem.get('idref') in ['Region.nonZeroRates', 'State.nonZeroRates']:
                to_remove.append(elem)
            if elem.tag == 'markovJumpsTreeLikelihood' and elem.get('idref') in ['Region.treeLikelihood', 'State.treeLikelihood']:
                to_remove.append(elem)
    
            # Spillover
            if elem.tag == 'column' and elem.get('label') == 'age(ingroup)':
                to_remove.append(elem)
            if elem.tag == 'column' and elem.get('label') == 'stem.proportion':
                to_remove.append(elem)
            if elem.tag == 'column' and elem.get('label') == 'age(transition)':
                to_remove.append(elem)
            if elem.tag == 'tmrcaStatistic' and elem.get('idref') in ['tmrca(ingroup)', 'tmrca(outgroup)', 'age(ingroup)', 'age(outgroup)']:
                to_remove.append(elem)
            if elem.tag == 'parameter' and elem.get('idref') == 'apobec3.stem.proportion':
                to_remove.append(elem)
            if elem.tag == 'parameter' and elem.get('idref') == 'apobec3.transition.time':
                to_remove.append(elem)
            if elem.tag == 'ageStatistic' and elem.get('idref') == 'age(apobec3.transition)':
                to_remove.append(elem)
            if elem.tag == 'parameter' and elem.get('idref') == 'constant.popSize':
                to_remove.append(elem)
            if (elem.tag == 'parameter' and
                elem.get('idref').startswith('age') and
                elem.get('idref') not in [f'age({s})' for s in included_sequences]):
                to_remove.append(elem)
            if elem.tag == 'ancestralTreeLikelihood' and elem.get('idref') == 'apobec3.treeLikelihood':
                elem.tag = 'treeDataLikelihood'
            if elem.tag == 'localClockModel' and elem.get('idref') == 'apobec3.branchRates':
                elem.tag = 'strictClockBranchRates'
            if elem.tag == 'coalescentLikelihood' and elem.get('idref') == 'coalescent.outgroup':
                to_remove.append(elem)
            if elem.tag == 'coalescentLikelihood' and elem.get('idref') == 'coalescent.ingroup':
                elem.set('idref', 'coalescent')
    
        for elem in to_remove:
            log.remove(elem)

    # Only modeling exponential growth post-spillover
    root.remove(root.find('constantSize'))
    startingTree = root.find('coalescentTree')
    assert startingTree.get('id') == 'startingTree'
    startingTree.find('taxa').set('idref', 'taxa')
    startingTree.find('constantSize').set('idref', 'exponential')
    startingTree.remove(startingTree.find('coalescentTree'))

    # And without DTA & spillover, we don't need very long chains at all to reach high ESSs
    mcmc.set('chainLength', '50000000') # was '400000000' in original file (8x longer!)

    beastXml.write(beast_path.as_posix())
