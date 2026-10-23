#!/usr/bin/env python3
import argparse
import gzip
import io
import logging
import shlex
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path
VERBOSE: bool = False
CURRENT_LOG_PATH: Path | None = None

def ui_step(num: int, total: int, title: str) -> None:
    print(f'[{num}/{total}] {title}')

def ui_info(msg: str) -> None:
    print(f'    {msg}')

def ui_ok(msg: str) -> None:
    print(f'    -> {msg}')

def ui_warn(msg: str) -> None:
    print(f'    [WARN] {msg}')

def ui_err(msg: str) -> None:
    print(f'    [ERROR] {msg}', file=sys.stderr)

def get_dataset_counts(prefix: str) -> tuple[int, int]:
    fam = Path(f'{prefix}.fam')
    bim = Path(f'{prefix}.bim')
    ind = Path(f'{prefix}.ind')
    snp = Path(f'{prefix}.snp')
    n_samples = 0
    n_variants = 0
    sample_file = fam if fam.is_file() else ind if ind.is_file() else None
    variant_file = bim if bim.is_file() else snp if snp.is_file() else None
    if sample_file and sample_file.is_file():
        try:
            with open(sample_file, 'r', encoding='utf-8', errors='ignore') as fh:
                n_samples = sum((1 for line in fh if line.strip()))
        except Exception:
            pass
    if variant_file and variant_file.is_file():
        try:
            with open(variant_file, 'r', encoding='utf-8', errors='ignore') as fh:
                n_variants = sum((1 for line in fh if line.strip()))
        except Exception:
            pass
    return (n_samples, n_variants)

def setup_logging(out_prefix: str, to_file: bool=True, verbose: bool=False) -> None:
    global VERBOSE, CURRENT_LOG_PATH
    VERBOSE = verbose
    handlers: list[logging.Handler] = []
    if verbose:
        handlers.append(logging.StreamHandler(sys.stdout))
    if to_file:
        log_path = Path(f'{out_prefix}_merge.log')
        log_path.parent.mkdir(parents=True, exist_ok=True)
        CURRENT_LOG_PATH = log_path
        handlers.append(logging.FileHandler(log_path, mode='w', encoding='utf-8'))
    if not handlers:
        handlers.append(logging.NullHandler())
    fmt = '%(asctime)s [%(levelname)s] %(message)s'
    logging.root.handlers = []
    logging.basicConfig(level=logging.DEBUG if to_file else logging.INFO, format=fmt, handlers=handlers)
    if to_file:
        logging.info(f'Log started: {out_prefix}_merge.log')

def check_tools(script_dir: Path) -> bool:
    print('=== GenoMerge Tool Dependency Check ===')
    all_ok = True
    tools = [('plink', 'PLINK 1.9', 'https://www.cog-genomics.org/plink/1.9/'), ('plink2', 'PLINK 2.0', 'https://www.cog-genomics.org/plink/2.0/')]
    for name, label, url in tools:
        found_path = None
        for candidate in [Path.cwd() / f'{name}.exe', Path.cwd() / name, script_dir / f'{name}.exe', script_dir / name]:
            if candidate.is_file():
                found_path = str(candidate.resolve())
                break
        if not found_path:
            found_path = shutil.which(name)
        if found_path:
            ver = 'detected'
            try:
                res = subprocess.run([found_path, '--version'], capture_output=True, text=True, timeout=5)
                lines = (res.stdout or res.stderr).strip().splitlines()
                if lines:
                    ver = lines[0]
            except Exception:
                pass
            print(f'  [OK] {label:10s} : {found_path}')
            print(f'       Version    : {ver}')
        else:
            all_ok = False
            print(f'  [MISSING] {label:10s} : Not found in current directory, package directory, or PATH')
            print(f'            Download  : {url}')
    if all_ok:
        print('\nAll required external tools are available.')
    else:
        print('\nPlease download the missing tool(s) and place them on your system PATH or in your working directory.')
    return all_ok

def find_exe(name: str, script_dir: Path) -> str:
    for candidate in [Path.cwd() / f'{name}.exe', Path.cwd() / name, script_dir / f'{name}.exe', script_dir / name]:
        if candidate.is_file():
            return str(candidate.resolve())
    found = shutil.which(name)
    if found:
        return found
    sys.exit(f"ERROR: '{name}' not found in current directory, {script_dir}, or PATH")

def run(cmd: list[str], quiet_fail: bool=False) -> int:
    cmd_str = ' '.join((str(c) for c in cmd))
    logging.info(f'RUN: {cmd_str}')
    result = subprocess.run(cmd, text=True, capture_output=True, encoding='utf-8', errors='replace')
    if result.stdout:
        logging.debug('STDOUT:\n' + result.stdout)
    if result.stderr:
        logging.debug('STDERR:\n' + result.stderr)
    if VERBOSE:
        if result.stdout:
            sys.stdout.write(result.stdout)
            sys.stdout.flush()
        if result.stderr:
            sys.stderr.write(result.stderr)
            sys.stderr.flush()
    if result.returncode != 0:
        logging.warning(f'Command exited {result.returncode}')
        if not VERBOSE and (not quiet_fail):
            combined = (result.stdout or '') + '\n' + (result.stderr or '')
            lines = [line.rstrip() for line in combined.splitlines() if line.strip()]
            tail = lines[-15:] if len(lines) > 15 else lines
            ui_err(f"Command exited {result.returncode}: {' '.join((str(c) for c in cmd[:4]))}...")
            for line in tail:
                print(f'        | {line}', file=sys.stderr)
            if CURRENT_LOG_PATH:
                print(f'        (Full details in: {CURRENT_LOG_PATH})', file=sys.stderr)
    return result.returncode

def open_raw_file(path: Path | str):
    p = Path(str(path).strip('"\''))
    if p.suffix == '.gz':
        return gzip.open(p, 'rt', encoding='utf-8', errors='ignore')
    elif p.suffix == '.zip':
        zf = zipfile.ZipFile(p)
        names = [n for n in zf.namelist() if not n.startswith('__') and n.lower().endswith(('.txt', '.csv', '.tsv', '.vcf'))]
        if not names:
            names = [n for n in zf.namelist() if not n.startswith('__')]
        if not names:
            raise FileNotFoundError(f'No usable text/genotype file found inside {p}')
        return io.TextIOWrapper(zf.open(names[0]), encoding='utf-8', errors='ignore')
    return open(p, 'r', encoding='utf-8', errors='ignore')

def is_vcf_file(path: Path) -> bool:
    p = Path(str(path).strip('"\''))
    if not p.is_file():
        return False
    name_l = p.name.lower()
    if name_l.endswith(('.vcf', '.vcf.gz')):
        return True
    try:
        with open_raw_file(p) as fh:
            for i, line in enumerate(fh):
                if i > 50:
                    break
                if line.startswith('##fileformat=VCF'):
                    return True
    except Exception:
        return False
    return False

def is_raw_dtc_file(path: Path) -> bool:
    p = Path(str(path).strip('"\''))
    if not p.is_file():
        return False
    chr_keys = set((str(i) for i in range(1, 27))) | {'X', 'Y', 'MT', 'M', 'XY'}
    try:
        with open_raw_file(p) as fh:
            for i, line in enumerate(fh):
                if i > 200:
                    break
                line_l = line.lower()
                if any((tag in line_l for tag in ('23andme', 'ancestrydna', 'myheritage', 'familytreedna', 'living dna'))):
                    return True
                if line.startswith('#') or not line.strip():
                    continue
                parts = [p.strip().strip('"\'') for p in (line.split(',') if ',' in line else line.split())]
                if not parts or parts[0].upper() in ('RSID', 'SNP', 'ID'):
                    continue
                chrom = parts[1].upper().removeprefix('CHR')
                if len(parts) in (4, 5) and parts[0].startswith(('rs', 'i', 'vg', 'kg')) and (chrom in chr_keys):
                    return True
    except Exception:
        return False
    return False

def strip_format_extension(path_str: str) -> str:
    p = Path(str(path_str).strip('"\''))
    stem = p.name
    for ext in ('.bed', '.bim', '.fam', '.geno', '.ind', '.snp'):
        if stem.lower().endswith(ext):
            clean_stem = stem[:-len(ext)]
            return str(p.parent / clean_stem) if str(p.parent) not in ('', '.') else clean_stem
    return str(path_str).strip('"\'')

def resolve_raw_path(prefix: str) -> Path | None:
    p = Path(str(prefix).strip('"\''))
    if p.is_file():
        return p
    for ext in ('.txt', '.csv', '.tsv', '.txt.gz', '.csv.gz', '.tsv.gz', '.zip', '.gz', '.vcf', '.vcf.gz'):
        cand = p.parent / f'{p.name}{ext}'
        if cand.is_file():
            return cand
    return None

def detect_format(prefix: str) -> str:
    p = Path(str(prefix).strip('"\''))
    parent, stem = (p.parent, p.name)
    if any((stem.lower().endswith(ext) for ext in ('.bed', '.bim', '.fam', '.geno', '.ind', '.snp'))):
        stem = stem.rsplit('.', 1)[0]
    if all(((parent / f'{stem}.{ext}').exists() for ext in ('bed', 'bim', 'fam'))):
        return 'plink'
    if all(((parent / f'{stem}.{ext}').exists() for ext in ('geno', 'ind', 'snp'))):
        return 'eigenstrat'
    resolved = resolve_raw_path(prefix)
    if resolved is not None:
        if is_vcf_file(resolved):
            return 'vcf'
        if is_raw_dtc_file(resolved):
            return 'raw'
    return 'error'

def clean_and_normalize_raw_dtc(input_path: Path, output_path: Path) -> str:
    chr_map = {str(i): i for i in range(1, 23)}
    chr_map.update({'X': 23, '23': 23, 'Y': 24, '24': 24, 'XY': 25, '25': 25, 'MT': 26, 'M': 26, '26': 26})
    provider = 'Commercial DTC Array'
    records = []
    with open_raw_file(input_path) as fh:
        for line in fh:
            line_str = line.strip()
            if not line_str:
                continue
            line_l = line_str.lower()
            if line_str.startswith('#'):
                if 'ancestrydna' in line_l:
                    provider = 'AncestryDNA'
                elif 'myheritage' in line_l:
                    provider = 'MyHeritage'
                elif '23andme' in line_l:
                    provider = '23andMe'
                elif 'living dna' in line_l:
                    provider = 'LivingDNA'
                elif 'familytreedna' in line_l:
                    provider = 'FamilyTreeDNA'
                continue
            if '\t' in line_str:
                parts = [p.strip().strip('"\'') for p in line_str.split('\t') if p.strip()]
            elif ',' in line_str:
                parts = [p.strip().strip('"\'') for p in line_str.split(',') if p.strip()]
            else:
                parts = [p.strip().strip('"\'') for p in line_str.split() if p.strip()]
            if not parts:
                continue
            if parts[0].upper() in ('RSID', 'SNP', 'ID'):
                if len(parts) >= 5:
                    provider = 'AncestryDNA / 5-col DTC'
                continue
            if len(parts) < 4:
                continue
            rsid = parts[0]
            chrom_str = parts[1].upper().removeprefix('CHR')
            pos_str = parts[2]
            if chrom_str not in chr_map:
                continue
            if len(parts) == 4:
                geno = parts[3].upper()
                if geno in ('00', '??', 'NN', '0', '?', '.', '-'):
                    geno = '--'
            else:
                a1, a2 = (parts[3].upper(), parts[4].upper())
                if a1 in ('0', '-', '?', '.', 'N') or a2 in ('0', '-', '?', '.', 'N'):
                    geno = '--'
                else:
                    geno = a1 + a2
            if geno == '--' or not any((c in 'ACGT' for c in geno)):
                continue
            try:
                pos = int(pos_str)
            except ValueError:
                continue
            records.append((chr_map[chrom_str], chrom_str, pos, rsid, geno))
    records.sort(key=lambda r: (r[0], r[2]))
    with open(output_path, 'w', encoding='utf-8') as out:
        out.write('# rsid\tchromosome\tposition\tgenotype\n')
        for _, chrom, pos, rsid, geno in records:
            out.write(f'{rsid}\t{chrom}\t{pos}\t{geno}\n')
    logging.info(f'Detected {provider} format: cleaned and sorted {len(records)} variants from {input_path}')
    return provider

def has_duplicate_variants(prefix: str) -> bool:
    bim_path = Path(f'{prefix}.bim')
    if not bim_path.is_file():
        return False
    seen = set()
    with open(bim_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                vid = parts[1]
                if vid != '.' and vid in seen:
                    return True
                seen.add(vid)
    return False

def convert_raw_dtc_to_plink(input_path: Path, out_prefix: str, fid: str, iid: str, sex: int, plink_exe: str, plink2_exe: str) -> str:
    clean_txt = Path(f'{out_prefix}_cleaned.txt')
    clean_and_normalize_raw_dtc(input_path, clean_txt)
    raw_plink = f'{out_prefix}_raw'
    rc = run([plink_exe, '--23file', str(clean_txt), fid, iid, str(sex), '0', '--snps-only', 'just-acgt', '--set-hh-missing', '--allow-no-sex', '--make-bed', '--out', raw_plink])
    if rc != 0 or not Path(f'{raw_plink}.bed').exists():
        sys.exit(f'ERROR: plink conversion failed for {input_path}')
    bim_path = Path(f'{raw_plink}.bim')
    bad_snps = []
    if bim_path.is_file():
        with open(bim_path, 'r', encoding='utf-8', errors='ignore') as fh:
            for line in fh:
                parts = line.strip().split()
                if len(parts) >= 6 and parts[4] == parts[5]:
                    bad_snps.append(parts[1])
    plink_for_dup = raw_plink
    bad_snps_file = None
    if bad_snps:
        bad_snps_file = Path(f'{out_prefix}_bad_alleles.snplist')
        bad_snps_file.write_text('\n'.join(bad_snps) + '\n', encoding='utf-8')
        clean_raw = f'{out_prefix}_raw_clean'
        run([plink_exe, '--bfile', raw_plink, '--exclude', str(bad_snps_file), '--set-hh-missing', '--allow-no-sex', '--make-bed', '--out', clean_raw], quiet_fail=True)
        if Path(f'{clean_raw}.bed').exists():
            plink_for_dup = clean_raw
    if has_duplicate_variants(plink_for_dup):
        logging.info(f'Duplicate variant IDs detected in {plink_for_dup} — deduplicating via plink2 --rm-dup...')
        final_plink = f'{out_prefix}_plink'
        rc = run([plink2_exe, '--bfile', plink_for_dup, '--rm-dup', 'force-first', '--make-bed', '--out', final_plink])
        if rc == 0 and Path(f'{final_plink}.bed').exists():
            if bad_snps_file and bad_snps_file.is_file():
                bad_snps_file.unlink(missing_ok=True)
            clean_txt.unlink(missing_ok=True)
            for ext in ('bed', 'bim', 'fam', 'log', 'nosex'):
                Path(f'{raw_plink}.{ext}').unlink(missing_ok=True)
                Path(f'{out_prefix}_raw_clean.{ext}').unlink(missing_ok=True)
            return final_plink
    if bad_snps_file and bad_snps_file.is_file():
        bad_snps_file.unlink(missing_ok=True)
    clean_txt.unlink(missing_ok=True)
    return plink_for_dup

def convert_vcf_to_plink(vcf_path: Path, out_prefix: str, fid: str, iid: str, sex: int, plink2_exe: str) -> str:
    final_plink = f'{out_prefix}_plink'
    rc = run([plink2_exe, '--vcf', str(vcf_path), '--double-id', '--snps-only', 'just-acgt', '--max-alleles', '2', '--make-bed', '--out', final_plink])
    if rc != 0 or not Path(f'{final_plink}.bed').exists():
        sys.exit(f'ERROR: plink2 VCF conversion failed for {vcf_path}')
    if has_duplicate_variants(final_plink):
        logging.info(f'Duplicate variant IDs detected in VCF {final_plink} — deduplicating via plink2 --rm-dup...')
        dedup_plink = f'{out_prefix}_dedup'
        rc = run([plink2_exe, '--bfile', final_plink, '--rm-dup', 'force-first', '--make-bed', '--out', dedup_plink])
        if rc == 0 and Path(f'{dedup_plink}.bed').exists():
            for ext in ('bed', 'bim', 'fam', 'log', 'nosex'):
                Path(f'{final_plink}.{ext}').unlink(missing_ok=True)
            final_plink = dedup_plink
    fam_path = Path(f'{final_plink}.fam')
    lines = [l for l in fam_path.read_text().splitlines() if l.strip()]
    is_single_sample = len(lines) == 1
    new_lines = []
    for line in lines:
        parts = line.split()
        if len(parts) >= 6:
            if fid:
                parts[0] = fid
            if iid and is_single_sample:
                parts[1] = iid
            if sex and is_single_sample:
                parts[4] = str(sex)
        new_lines.append('\t'.join(parts))
    fam_path.write_text('\n'.join(new_lines) + '\n')
    return final_plink

def bk_col(prefix: str) -> None:
    shutil.copy2(f'{prefix}.ind', f'{prefix}_backup.ind')
    logging.info(f'Backed up {prefix}.ind -> {prefix}_backup.ind')

def nb_col(prefix: str) -> None:
    lines = [l for l in Path(f'{prefix}.ind').read_text().splitlines() if l.strip()]
    new_lines = []
    for i, line in enumerate(lines, 1):
        parts = line.split()
        if len(parts) >= 3:
            parts[0] = str(i)
            parts[2] = str(i)
        new_lines.append('\t'.join(parts))
    Path(f'{prefix}.ind').write_text('\n'.join(new_lines) + '\n')
    logging.info(f'Numbered col1/col3 in {prefix}.ind')

def ps_col(plink_prefix: str, orig_prefix: str) -> None:
    backup_lines = [l for l in Path(f'{orig_prefix}_backup.ind').read_text().splitlines() if l.strip()]
    fam_path = Path(f'{plink_prefix}.fam')
    fam_lines = [l for l in fam_path.read_text().splitlines() if l.strip()]
    if len(backup_lines) != len(fam_lines):
        sys.exit(f'ERROR: sample count mismatch — {orig_prefix}_backup.ind ({len(backup_lines)}) vs {plink_prefix}.fam ({len(fam_lines)})')
    new_fam = []
    for ind_line, fam_line in zip(backup_lines, fam_lines):
        ind_parts = ind_line.split()
        fam_parts = fam_line.split()
        fam_parts[0] = ind_parts[2]
        fam_parts[1] = ind_parts[0]
        if len(fam_parts) >= 6:
            fam_parts[5] = '-9'
        new_fam.append('\t'.join(fam_parts))
    fam_path.write_text('\n'.join(new_fam) + '\n')
    logging.info(f'Restored original IDs into {plink_prefix}.fam')

def convert_to_plink(prefix: str, plink2_exe: str, out_prefix: str | None=None) -> str:
    out = out_prefix or f'{prefix}_plink'
    rc = run([plink2_exe, '--eigfile', prefix, 'nohash', '--make-bed', '--out', out])
    if rc != 0 or not Path(f'{out}.bed').exists():
        sys.exit(f'ERROR: plink2 conversion failed for {prefix}')
    return out

def plink_merge(ref: str, sub: str, out: str, plink_exe: str) -> int:
    Path(f'{out}-merge.missnp').unlink(missing_ok=True)
    return run([plink_exe, '--bfile', ref, '--bmerge', sub, '--keep-allele-order', '--indiv-sort', '0', '--set-hh-missing', '--allow-no-sex', '--make-bed', '--out', out], quiet_fail=True)

def has_missnp(out: str) -> bool:
    return Path(f'{out}-merge.missnp').exists()

def merge_pipeline(ref: str, sub: str, out: str, plink_exe: str, step_prefix: str | None=None, plink2_exe: str | None=None, chr_filter: str | None=None) -> bool:
    pfx = step_prefix or out
    sub_vars = get_bim_variants(sub)
    ref_vars = get_bim_variants(ref)
    if sub_vars and sub_vars.issubset(ref_vars):
        logging.info(f"Subject '{sub}' verified as strict subset of reference ({len(sub_vars):,} variants)")
        sub_target = sub
    else:
        snplist = f'{pfx}_Ref_snp'
        snplist_cmd = [plink_exe, '--bfile', ref]
        if chr_filter:
            snplist_cmd.extend(['--chr', chr_filter])
        snplist_cmd.extend(['--snps-only', 'just-acgt', '--write-snplist', '--allow-no-sex', '--out', snplist])
        rc = run(snplist_cmd)
        if rc != 0:
            logging.error(f'--write-snplist failed for {ref}')
            return False
        sub_ext = f'{pfx}_Sub_extracted'
        sub_ext_cmd = [plink_exe, '--bfile', sub, '--extract', f'{snplist}.snplist']
        if chr_filter:
            sub_ext_cmd.extend(['--chr', chr_filter])
        sub_ext_cmd.extend(['--snps-only', 'just-acgt', '--set-hh-missing', '--allow-no-sex', '--make-bed', '--out', sub_ext])
        rc = run(sub_ext_cmd)
        if rc != 0 or not Path(f'{sub_ext}.bed').exists():
            logging.error(f'--extract failed for {sub}')
            return False
        sub_target = sub_ext
    if plink2_exe and has_duplicate_variants(sub_target):
        logging.info(f'Duplicate variant IDs detected in subject set — deduplicating via plink2 --rm-dup...')
        sub_dedup = f'{pfx}_Sub_dedup'
        rc = run([plink2_exe, '--bfile', sub_target, '--rm-dup', 'force-first', '--make-bed', '--out', sub_dedup], quiet_fail=True)
        if rc == 0 and Path(f'{sub_dedup}.bed').exists():
            sub_target = sub_dedup
    logging.info('Merge attempt 1 (straight merge)...')
    plink_merge(ref, sub_target, out, plink_exe)
    if not has_missnp(out) and Path(f'{out}.bed').exists():
        ns, nv = get_dataset_counts(out)
        ui_ok(f'Merged successfully on attempt 1 ({ns:,} samples, {nv:,} variants)')
        logging.info('SUCCESS: merge completed (attempt 1)')
        return True
    ui_info('Strand conflicts detected -> Attempt 2: flipping missnp variants...')
    logging.info('Merge attempt 2 (flip missnp)...')
    sub_flip = f'{pfx}_Sub_flipped'
    run([plink_exe, '--bfile', sub_target, '--flip', f'{out}-merge.missnp', '--set-hh-missing', '--allow-no-sex', '--make-bed', '--out', sub_flip])
    sub_attempt2 = sub_flip if Path(f'{sub_flip}.bed').exists() else sub_target
    plink_merge(ref, sub_attempt2, out, plink_exe)
    if not has_missnp(out) and Path(f'{out}.bed').exists():
        ns, nv = get_dataset_counts(out)
        ui_ok(f'Merged successfully on attempt 2 after flip ({ns:,} samples, {nv:,} variants)')
        logging.info('SUCCESS: merge completed (attempt 2, after flip)')
        return True
    ui_info('Residual conflicts detected -> Attempt 3: excluding unresolvable missnp variants...')
    logging.info('Merge attempt 3 (exclude missnp)...')
    sub_clean = f'{pfx}_Sub_clean'
    sub_to_clean = sub_flip if Path(f'{sub_flip}.bed').exists() else sub_target
    run([plink_exe, '--bfile', sub_to_clean, '--exclude', f'{out}-merge.missnp', '--set-hh-missing', '--allow-no-sex', '--make-bed', '--out', sub_clean])
    sub_attempt3 = sub_clean if Path(f'{sub_clean}.bed').exists() else sub_to_clean
    plink_merge(ref, sub_attempt3, out, plink_exe)
    if not has_missnp(out) and Path(f'{out}.bed').exists():
        ns, nv = get_dataset_counts(out)
        ui_ok(f'Merged successfully on attempt 3 after exclude ({ns:,} samples, {nv:,} variants)')
        logging.info('SUCCESS: merge completed (attempt 3, after exclude)')
        return True
    ui_err(f'Merge failed after 3 attempts for {sub} into {ref}')
    logging.error(f'ERROR: merge failed after 3 attempts for {sub} into {ref}')
    return False

def export_eigenstrat(out: str, plink2_exe: str) -> None:
    logging.info('Exporting merged dataset to eigenstrat format natively (preserves ADMIXTOOLS hash)...')
    fam_path = Path(f'{out}.fam')
    if not fam_path.is_file():
        sys.exit(f'ERROR: {fam_path} not found for EIGENSTRAT export')
    fam_lines = [l for l in fam_path.read_text(encoding='utf-8', errors='ignore').splitlines() if l.strip()]
    sample_metadata = []
    clean_fam = []
    for line in fam_lines:
        parts = line.split()
        if len(parts) >= 2:
            fid = parts[0]
            iid = parts[1]
            pid = parts[2] if len(parts) > 2 else '0'
            mid = parts[3] if len(parts) > 3 else '0'
            raw_sex = parts[4] if len(parts) > 4 else '0'
            sex = 'M' if raw_sex in ('1', 'M', 'm') else ('F' if raw_sex in ('2', 'F', 'f') else 'U')
            clean_iid = iid[:30] + '~' + iid[-8:] if len(iid) > 39 else iid
            clean_pop = fid[:30] + '~' + fid[-8:] if len(fid) > 39 else fid
            sample_metadata.append((clean_iid, sex, clean_pop))
            clean_fam.append(f"0\t{clean_iid}\t{pid}\t{mid}\t{raw_sex if raw_sex in ('1', '2') else '0'}\t-9")
    fam_path.write_text('\n'.join(clean_fam) + '\n', encoding='utf-8')
    rc = run([plink2_exe, '--bfile', out, '--no-categorical', '--snps-only', 'just-acgt', '--export', 'eig', 'id-paste=iid', '--out', out])
    if rc != 0 or not Path(f'{out}.geno').exists():
        sys.exit(f'ERROR: eigenstrat export failed for {out}')
    ind_lines = [f'{s_iid}\t{s_sex}\t{s_pop}' for s_iid, s_sex, s_pop in sample_metadata]
    Path(f'{out}.ind').write_text('\n'.join(ind_lines) + '\n', encoding='utf-8')
    logging.info(f'Generated {out}.geno/.snp/.ind natively with valid EIGENSTRAT hash.')

def format_plink_fam(fam_path: Path) -> None:
    if not fam_path.is_file():
        return
    lines = [l for l in fam_path.read_text(encoding='utf-8', errors='ignore').splitlines() if l.strip()]
    cleaned = []
    for line in lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        fid = parts[0]
        iid = parts[1]
        pid = parts[2] if len(parts) > 2 and parts[2] not in ('.', 'NA') else '0'
        mid = parts[3] if len(parts) > 3 and parts[3] not in ('.', 'NA') else '0'
        sex = '0'
        if len(parts) > 4:
            s = parts[4].upper()
            if s in ('1', 'M', 'MALE'):
                sex = '1'
            elif s in ('2', 'F', 'FEMALE'):
                sex = '2'
            elif s in ('0', 'U', 'UNKNOWN', '-9'):
                sex = '0'
            elif s.isdigit():
                sex = s
        pheno = '-9'
        if len(parts) > 5 and parts[5] not in ('0', '-9', '-', 'NA', '.'):
            pheno = parts[5]
        cleaned.append(f'{fid}\t{iid}\t{pid}\t{mid}\t{sex}\t{pheno}')
    fam_path.write_text('\n'.join(cleaned) + '\n', encoding='utf-8')
    logging.info(f'Standardized {fam_path} to 6-column PLINK format ({len(cleaned)} samples).')

def cleanup_step(step_prefix: str, step_out: str, prev_ref: str | None=None, temp_sub: str | None=None) -> None:
    p_step = Path(step_prefix)
    parent = p_step.parent if str(p_step.parent) != '' else Path('.')
    stem = p_step.name
    pats = [f'{stem}_Ref_snp*', f'{stem}_Sub_extracted*', f'{stem}_Sub_dedup*', f'{stem}_Sub_flipped*', f'{stem}_Sub_clean*']
    p_out = Path(step_out)
    out_parent = p_out.parent if str(p_out.parent) != '' else Path('.')
    out_stem = p_out.name
    pats_out = [f'{out_stem}-merge*']
    for pat in pats:
        for f in parent.glob(pat):
            f.unlink(missing_ok=True)
            logging.info(f'Removed step temp: {f}')
    for pat in pats_out:
        for f in out_parent.glob(pat):
            f.unlink(missing_ok=True)
            logging.info(f'Removed step temp: {f}')
    if prev_ref is not None:
        p_prev = Path(prev_ref)
        prev_parent = p_prev.parent if str(p_prev.parent) != '' else Path('.')
        prev_stem = p_prev.name
        for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
            target = prev_parent / f'{prev_stem}{ext}'
            if target.is_file():
                target.unlink(missing_ok=True)
                logging.info(f'Removed prev step ref: {target}')
    if temp_sub is not None:
        p_sub = Path(temp_sub)
        sub_parent = p_sub.parent if str(p_sub.parent) != '' else Path('.')
        sub_stem = p_sub.name
        for ext in ('.bed', '.bim', '.fam', '.log', '.nosex', '.txt'):
            target = sub_parent / f'{sub_stem}{ext}'
            if target.is_file():
                target.unlink(missing_ok=True)
                logging.info(f'Removed temp subject: {target}')

def cleanup(out: str, delete_final_plink: bool=False, extra_prefixes: list[str] | None=None) -> None:
    out_p = Path(out)
    parent_dir = out_p.parent if str(out_p.parent) != '' else Path('.')
    out_stem = out_p.name
    patterns = [f'{out_stem}_Ref_snp*', f'{out_stem}_ref_snplist*', f'{out_stem}_ref_filtered*', f'{out_stem}_Sub_*', f'{out_stem}-merge*', f'{out_stem}_subraw*', f'{out_stem}_subvcf*', f'{out_stem}_sub23*', f'{out_stem}_step*', f'{out_stem}_ref_plink*', f'{out_stem}_norm_*']
    if delete_final_plink:
        patterns.extend([f'{out_stem}.bed', f'{out_stem}.bim', f'{out_stem}.fam', f'{out_stem}.nosex', f'{out_stem}.log'])
    for pattern in patterns:
        for f in parent_dir.glob(pattern):
            f.unlink(missing_ok=True)
            logging.info(f'Removed temp: {f}')
    if extra_prefixes:
        for ep in extra_prefixes:
            p_ep = Path(ep)
            ep_parent = p_ep.parent if str(p_ep.parent) != '' else Path('.')
            ep_stem = p_ep.name
            for ext in ('.bed', '.bim', '.fam', '.log', '.nosex', '.snplist', '-merge.missnp', '.txt'):
                target = ep_parent / f'{ep_stem}{ext}'
                if target.is_file():
                    target.unlink(missing_ok=True)
                    logging.info(f'Removed extra temp: {target}')

def get_dataset_stem(path_str: str) -> str:
    p = Path(path_str)
    name = p.name
    for sfx in ('.vcf.gz', '.txt.gz', '.csv.gz', '.tsv.gz'):
        if name.lower().endswith(sfx):
            return name[:-len(sfx)]
    return p.stem

def parse_batch_file(batch_path: Path) -> list[tuple[str, str | None]]:
    if not batch_path.is_file():
        sys.exit(f"ERROR: Batch file not found: '{batch_path}'")
    known_formats = {'plink', 'bed', 'eigenstrat', 'eig', 'vcf', 'raw', 'dtc', '23andme', '23'}
    entries: list[tuple[str, str | None]] = []
    with open(batch_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            line_str = line.strip()
            if not line_str or line_str.startswith('#'):
                continue
            parts = [p.strip().strip('"\'') for p in (line_str.split(',') if ',' in line_str else line_str.split()) if p.strip()]
            if not parts:
                continue
            if len(parts) == 1:
                entries.append((parts[0], None))
            elif len(parts) >= 2:
                if parts[-1].lower() in known_formats:
                    fmt_hint = parts[-1].lower()
                    ds_path = ' '.join(parts[:-1]) if ',' not in line_str else parts[0]
                    entries.append((ds_path, fmt_hint))
                elif parts[0].lower() in known_formats:
                    fmt_hint = parts[0].lower()
                    ds_path = ' '.join(parts[1:]) if ',' not in line_str else parts[1]
                    entries.append((ds_path, fmt_hint))
                else:
                    entries.append((line_str, None))
    return entries

def get_fam_samples(fam_path: Path) -> list[tuple[str, str]]:
    samples = []
    if not fam_path.is_file():
        return samples
    with open(fam_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                samples.append((parts[0], parts[1]))
    return samples

def check_duplicate_samples(datasets: list[tuple[str, str]]) -> list[tuple[tuple[str, str], list[str]]]:
    seen: dict[tuple[str, str], list[str]] = {}
    for label, prefix in datasets:
        fam_path = Path(f'{prefix}.fam')
        for sample in get_fam_samples(fam_path):
            seen.setdefault(sample, []).append(label)
    duplicates = [(sample, labels) for sample, labels in seen.items() if len(labels) > 1]
    return duplicates

def get_bim_variants(prefix: str) -> set[str]:
    bim_path = Path(f'{prefix}.bim')
    variants = set()
    if not bim_path.is_file():
        return variants
    with open(bim_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                variants.add(parts[1])
    return variants
AUTOSOMAL_CHRS = {str(i) for i in range(1, 23)}

def bim_has_other_chromosomes(prefix: str, allowed_chrs: set[str]) -> bool:
    bim_path = Path(f'{prefix}.bim')
    if not bim_path.is_file():
        return True
    with open(bim_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            parts = line.strip().split(maxsplit=1)
            if parts and parts[0] not in allowed_chrs:
                return True
    return False

def load_id_list(inputs: list[str] | str | None) -> set[str]:
    if not inputs:
        return set()
    if isinstance(inputs, str):
        inputs = [inputs]
    ids = set()
    for item in inputs:
        p = Path(item)
        if p.is_file():
            with open(p, 'r', encoding='utf-8', errors='ignore') as fh:
                for line in fh:
                    line_str = line.strip()
                    if line_str and (not line_str.startswith('#')):
                        parts = line_str.split()
                        if parts:
                            ids.add(parts[0])
        else:
            item_lower = item.lower()
            file_extensions = ('.txt', '.poplist', '.list', '.tsv', '.csv', '.fam', '.ind')
            if any((c in item for c in ('/', '\\'))) or any((item_lower.endswith(ext) for ext in file_extensions)):
                sys.exit(f"ERROR: file not found: '{item}'")
            if item.strip() and (not item.startswith('#')):
                ids.add(item.strip())
    return ids

def filter_plink_samples(bfile_prefix: str, out_prefix: str, plink_exe: str, target_fids: set[str], target_iids: set[str], dataset_name: str='reference', allow_empty: bool=False) -> tuple[str | None, bool]:
    fam_path = Path(f'{bfile_prefix}.fam')
    if not fam_path.is_file():
        sys.exit(f"ERROR: Cannot find .fam file for '{bfile_prefix}'")
    matching_samples: list[tuple[str, str]] = []
    total_samples = 0
    with open(fam_path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            parts = line.strip().split()
            if len(parts) >= 2:
                total_samples += 1
                fid, iid = (parts[0], parts[1])
                if target_fids and fid not in target_fids:
                    continue
                if target_iids and iid not in target_iids:
                    continue
                matching_samples.append((fid, iid))
    if not matching_samples:
        if allow_empty:
            logging.warning(f'0 of {total_samples} samples in {dataset_name} matched keep criteria.')
            return (None, False)
        sys.exit(f"ERROR: 0 of {total_samples} samples in {dataset_name} ('{bfile_prefix}') matched the keep criteria.")
    if len(matching_samples) == total_samples:
        logging.info(f'All {total_samples} samples in {dataset_name} matched keep criteria — keeping all.')
        return (bfile_prefix, False)
    keep_file = Path(f'{out_prefix}_keep.txt')
    with open(keep_file, 'w', encoding='utf-8') as out_fh:
        for fid, iid in matching_samples:
            out_fh.write(f'{fid}\t{iid}\n')
    filtered_plink = f'{out_prefix}'
    rc = run([plink_exe, '--bfile', bfile_prefix, '--keep', str(keep_file), '--make-bed', '--allow-no-sex', '--out', filtered_plink])
    keep_file.unlink(missing_ok=True)
    if rc != 0 or not Path(f'{filtered_plink}.bed').exists():
        sys.exit(f'ERROR: plink sample filtering failed for {dataset_name} ({bfile_prefix})')
    logging.info(f'Filtered {dataset_name} from {total_samples} to {len(matching_samples)} samples matching keep criteria.')
    ui_ok(f'Filtered {dataset_name}: kept {len(matching_samples):,} of {total_samples:,} samples')
    return (filtered_plink, True)

def check_king_cutoff(bfile_prefix: str, cutoff: float, plink2_exe: str, temp_prefix: str) -> list[tuple[str, str]]:
    king_out = f'{temp_prefix}_king_chk'
    rc = run([plink2_exe, '--bfile', bfile_prefix, '--king-cutoff', str(cutoff), '--out', king_out], quiet_fail=True)
    removed: list[tuple[str, str]] = []
    out_id_file = Path(f'{king_out}.king.cutoff.out.id')
    if out_id_file.is_file():
        with open(out_id_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line_str = line.strip()
                if line_str and (not line_str.startswith('#')):
                    parts = line_str.split()
                    if len(parts) >= 2:
                        removed.append((parts[0], parts[1]))
    for ext in ('.log', '.king.cutoff.out.id', '.king.cutoff.in.id'):
        Path(f'{king_out}{ext}').unlink(missing_ok=True)
    return removed

def apply_king_cutoff(bfile_prefix: str, out_prefix: str, cutoff: float, plink2_exe: str, plink_exe: str, dataset_label: str='dataset') -> tuple[str, list[tuple[str, str]], bool]:
    king_out = f'{out_prefix}_king'
    rc = run([plink2_exe, '--bfile', bfile_prefix, '--king-cutoff', str(cutoff), '--out', king_out], quiet_fail=True)
    out_id_file = Path(f'{king_out}.king.cutoff.out.id')
    in_id_file = Path(f'{king_out}.king.cutoff.in.id')
    if rc != 0 or not in_id_file.is_file():
        logging.warning(f'plink2 --king-cutoff failed or produced no output for {dataset_label}.')
        return (bfile_prefix, [], False)
    removed: list[tuple[str, str]] = []
    if out_id_file.is_file():
        with open(out_id_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                line_str = line.strip()
                if line_str and (not line_str.startswith('#')):
                    parts = line_str.split()
                    if len(parts) >= 2:
                        removed.append((parts[0], parts[1]))
    Path(f'{king_out}.log').unlink(missing_ok=True)
    Path(f'{king_out}.king.cutoff.out.id').unlink(missing_ok=True)
    if not removed:
        Path(f'{king_out}.king.cutoff.in.id').unlink(missing_ok=True)
        return (bfile_prefix, [], False)
    filtered_pfx = f'{out_prefix}_king_filtered'
    rc_keep = run([plink_exe, '--bfile', bfile_prefix, '--keep', str(in_id_file), '--allow-no-sex', '--make-bed', '--out', filtered_pfx], quiet_fail=True)
    Path(f'{king_out}.king.cutoff.in.id').unlink(missing_ok=True)
    if rc_keep != 0 or not Path(f'{filtered_pfx}.bed').exists():
        logging.warning(f'PLINK --keep failed after king-cutoff for {dataset_label}. Proceeding with unfiltered.')
        return (bfile_prefix, [], False)
    return (filtered_pfx, removed, True)

def apply_qc_filters(bfile_prefix: str, out_prefix: str, plink_exe: str, geno: float | None=None, mind: float | None=None, maf: float | None=None, mac: int | None=None, exclude_high_ld: bool=False) -> tuple[str, bool]:
    has_any = geno is not None or mind is not None or maf is not None or (mac is not None) or exclude_high_ld
    if not has_any:
        return (bfile_prefix, False)
    curr_pfx = bfile_prefix
    is_temp = False
    if exclude_high_ld:
        range_file = Path(f'{out_prefix}_high_ld_ranges.txt')
        lines = ['6\t25000000\t35000000\tHLA', '8\t7000000\t13000000\t8p23_inversion', '17\t40000000\t45000000\t17q21_inversion', '2\t135000000\t138000000\tLCT']
        range_file.write_text('\n'.join(lines) + '\n', encoding='utf-8')
        hld_out = f'{out_prefix}_no_high_ld'
        rc = run([plink_exe, '--bfile', curr_pfx, '--exclude', 'range', str(range_file), '--allow-no-sex', '--make-bed', '--out', hld_out], quiet_fail=True)
        range_file.unlink(missing_ok=True)
        if rc == 0 and Path(f'{hld_out}.bed').exists():
            if is_temp:
                for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                    Path(f'{curr_pfx}{ext}').unlink(missing_ok=True)
            curr_pfx = hld_out
            is_temp = True
            n_s, n_v = get_dataset_counts(curr_pfx)
            ui_ok(f'Excluded high-LD inversion regions (HLA, 8p23, 17q21, LCT) -> {n_v:,} variants remaining')
        else:
            logging.warning('Failed to exclude high-LD regions; continuing without it.')
    qc_args = []
    labels = []
    if mind is not None:
        qc_args.extend(['--mind', str(mind)])
        labels.append(f'mind={mind}')
    if geno is not None:
        qc_args.extend(['--geno', str(geno)])
        labels.append(f'geno={geno}')
    if maf is not None:
        qc_args.extend(['--maf', str(maf)])
        labels.append(f'maf={maf}')
    if mac is not None:
        qc_args.extend(['--mac', str(mac)])
        labels.append(f'mac={mac}')
    if qc_args:
        qc_out = f'{out_prefix}_qc'
        cmd = [plink_exe, '--bfile', curr_pfx, '--allow-no-sex'] + qc_args + ['--make-bed', '--out', qc_out]
        rc = run(cmd, quiet_fail=True)
        if rc == 0 and Path(f'{qc_out}.bed').exists():
            if is_temp:
                for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                    Path(f'{curr_pfx}{ext}').unlink(missing_ok=True)
            curr_pfx = qc_out
            is_temp = True
            n_s, n_v = get_dataset_counts(curr_pfx)
            ui_ok(f"QC filters applied ({', '.join(labels)}) -> {n_s:,} samples, {n_v:,} variants")
        else:
            logging.warning('QC filtering step failed; continuing with previous dataset.')
    return (curr_pfx, is_temp)

def apply_ld_pruning(bfile_prefix: str, out_prefix: str, plink_exe: str, window: int, step: int, r2: float) -> tuple[str, bool]:
    prune_tmp = f'{out_prefix}_prune_calc'
    rc1 = run([plink_exe, '--bfile', bfile_prefix, '--indep-pairwise', str(window), str(step), str(r2), '--allow-no-sex', '--out', prune_tmp], quiet_fail=True)
    in_file = Path(f'{prune_tmp}.prune.in')
    if rc1 != 0 or not in_file.is_file():
        logging.warning('LD pruning calculation failed; proceeding with unpruned data.')
        return (bfile_prefix, False)
    n_pruned = 0
    try:
        with open(in_file, 'r', encoding='utf-8', errors='ignore') as fh:
            n_pruned = sum((1 for line in fh if line.strip()))
    except Exception:
        pass
    ld_out = f'{out_prefix}_ld_pruned'
    rc2 = run([plink_exe, '--bfile', bfile_prefix, '--extract', str(in_file), '--allow-no-sex', '--make-bed', '--out', ld_out], quiet_fail=True)
    for ext in ('.prune.in', '.prune.out', '.log', '.nosex'):
        Path(f'{prune_tmp}{ext}').unlink(missing_ok=True)
    if rc2 == 0 and Path(f'{ld_out}.bed').exists():
        ui_ok(f'LD pruning ({window} {step} {r2}): retained {n_pruned:,} independent SNPs')
        return (ld_out, True)
    else:
        logging.warning('Failed to extract LD pruned SNPs; proceeding with unpruned data.')
        return (bfile_prefix, False)

def apply_variant_filters(bfile_prefix: str, out_prefix: str, plink_exe: str, extract_file: str | None=None, exclude_file: str | None=None) -> tuple[str, bool]:
    if not extract_file and not exclude_file:
        return (bfile_prefix, False)
    curr_pfx = bfile_prefix
    is_temp = False
    if extract_file:
        p_ext = Path(extract_file)
        if not p_ext.is_file():
            sys.exit(f"ERROR: --extract file '{extract_file}' not found")
        ext_out = f'{out_prefix}_extracted'
        rc = run([plink_exe, '--bfile', curr_pfx, '--extract', str(p_ext), '--allow-no-sex', '--make-bed', '--out', ext_out], quiet_fail=True)
        if rc == 0 and Path(f'{ext_out}.bed').exists():
            if is_temp:
                for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                    Path(f'{curr_pfx}{ext}').unlink(missing_ok=True)
            curr_pfx = ext_out
            is_temp = True
            n_s, n_v = get_dataset_counts(curr_pfx)
            ui_ok(f"Extracted variants matching '{p_ext.name}' -> {n_v:,} variants remaining")
            logging.info(f"Applied --extract with '{extract_file}': {n_v:,} variants retained")
        else:
            logging.warning(f"Failed to apply --extract with '{extract_file}'; proceeding without it.")
    if exclude_file:
        p_exc = Path(exclude_file)
        if not p_exc.is_file():
            sys.exit(f"ERROR: --exclude file '{exclude_file}' not found")
        exc_out = f'{out_prefix}_excluded'
        rc = run([plink_exe, '--bfile', curr_pfx, '--exclude', str(p_exc), '--allow-no-sex', '--make-bed', '--out', exc_out], quiet_fail=True)
        if rc == 0 and Path(f'{exc_out}.bed').exists():
            if is_temp:
                for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                    Path(f'{curr_pfx}{ext}').unlink(missing_ok=True)
            curr_pfx = exc_out
            is_temp = True
            n_s, n_v = get_dataset_counts(curr_pfx)
            ui_ok(f"Excluded variants matching '{p_exc.name}' -> {n_v:,} variants remaining")
            logging.info(f"Applied --exclude with '{exclude_file}': {n_v:,} variants retained")
        else:
            logging.warning(f"Failed to apply --exclude with '{exclude_file}'; proceeding without it.")
    return (curr_pfx, is_temp)

def normalize_ref_to_plink(ref_raw: str, out_prefix: str, plink2_exe: str) -> tuple[str, bool]:
    fmt_a = detect_format(ref_raw)
    if fmt_a == 'error':
        sys.exit(f"ERROR: cannot find PLINK or eigenstrat files for prefix '{ref_raw}'")
    if fmt_a == 'eigenstrat':
        logging.info(f'{ref_raw} is eigenstrat — converting to PLINK binary')
        bk_col(ref_raw)
        nb_col(ref_raw)
        ref_plink = convert_to_plink(ref_raw, plink2_exe, out_prefix=f'{out_prefix}_ref_plink')
        ps_col(ref_plink, ref_raw)
        shutil.copy2(f'{ref_raw}_backup.ind', f'{ref_raw}.ind')
        Path(f'{ref_raw}_backup.ind').unlink(missing_ok=True)
        return (ref_plink, True)
    elif fmt_a == 'plink':
        p = Path(ref_raw)
        if p.name.endswith(('.bed', '.bim', '.fam')):
            plink_stem = str(p.parent / p.stem)
        else:
            plink_stem = ref_raw
        logging.info(f'{ref_raw} is PLINK binary — using as-is')
        return (plink_stem, False)
    else:
        sys.exit(f"ERROR: reference '{ref_raw}' must be PLINK binary or EIGENSTRAT")

def normalize_subject_to_plink(sub_raw: str, format_hint: str | None, out_prefix: str, default_fid: str | None, default_iid: str | None, sex: int, plink_exe: str, plink2_exe: str, ref_snplist_path: Path | None=None, chr_filter: str | None=None) -> tuple[str, bool]:
    resolved_path = resolve_raw_path(sub_raw) or Path(sub_raw)
    fmt = None
    if format_hint:
        h = format_hint.lower()
        if h in ('vcf',):
            fmt = 'vcf'
        elif h in ('raw', 'dtc', '23andme', '23'):
            fmt = 'raw'
        elif h in ('eigenstrat', 'eig'):
            fmt = 'eigenstrat'
        elif h in ('plink', 'bed'):
            fmt = 'plink'
    if fmt is None:
        fmt = detect_format(sub_raw)
    stem_name = get_dataset_stem(str(resolved_path))
    sub_fid = default_fid or stem_name
    sub_iid = default_iid or stem_name
    converted = None
    is_temp = False
    if fmt == 'vcf':
        if not resolved_path.is_file():
            sys.exit(f"ERROR: VCF file not found: '{sub_raw}'")
        logging.info(f'{resolved_path} is VCF format — converting to PLINK binary')
        converted = convert_vcf_to_plink(resolved_path, out_prefix, sub_fid, sub_iid, sex, plink2_exe)
        is_temp = True
    elif fmt in ('raw', '23andme'):
        if not resolved_path.is_file():
            sys.exit(f"ERROR: raw DTC file not found: '{sub_raw}'")
        logging.info(f'{resolved_path} is raw commercial array format — normalizing and converting to PLINK binary')
        converted = convert_raw_dtc_to_plink(resolved_path, out_prefix, sub_fid, sub_iid, sex, plink_exe, plink2_exe)
        is_temp = True
    elif fmt == 'eigenstrat':
        logging.info(f'{sub_raw} is eigenstrat — converting to PLINK binary')
        bk_col(sub_raw)
        nb_col(sub_raw)
        converted = convert_to_plink(sub_raw, plink2_exe, out_prefix=out_prefix)
        ps_col(converted, sub_raw)
        shutil.copy2(f'{sub_raw}_backup.ind', f'{sub_raw}.ind')
        Path(f'{sub_raw}_backup.ind').unlink(missing_ok=True)
        is_temp = True
    elif fmt == 'plink':
        p = Path(sub_raw)
        if p.name.endswith(('.bed', '.bim', '.fam')):
            converted = str(p.parent / p.stem)
        else:
            converted = sub_raw
        logging.info(f'{sub_raw} is PLINK binary')
        is_temp = False
    else:
        sys.exit(f"ERROR: cannot find PLINK, eigenstrat, raw DTC, or VCF files for '{sub_raw}'")
    if ref_snplist_path and ref_snplist_path.is_file():
        logging.info(f"Extracting strictly reference SNP subset for '{stem_name}'...")
        subset_prefix = f'{out_prefix}_subset'
        subset_cmd = [plink_exe, '--bfile', converted, '--extract', str(ref_snplist_path)]
        if chr_filter:
            subset_cmd.extend(['--chr', chr_filter])
        subset_cmd.extend(['--snps-only', 'just-acgt', '--set-hh-missing', '--allow-no-sex', '--make-bed', '--out', subset_prefix])
        rc = run(subset_cmd, quiet_fail=True)
        if rc != 0 or not Path(f'{subset_prefix}.bed').exists():
            return (converted, is_temp)
        if is_temp:
            for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                Path(f'{converted}{ext}').unlink(missing_ok=True)
        final_pfx = subset_prefix
        if has_duplicate_variants(subset_prefix):
            logging.info(f'Duplicate variant IDs detected in subset {subset_prefix} — deduplicating via plink2 --rm-dup...')
            dedup_prefix = f'{out_prefix}_dedup'
            rc_dup = run([plink2_exe, '--bfile', subset_prefix, '--rm-dup', 'force-first', '--make-bed', '--out', dedup_prefix], quiet_fail=True)
            if rc_dup == 0 and Path(f'{dedup_prefix}.bed').exists():
                for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                    Path(f'{subset_prefix}{ext}').unlink(missing_ok=True)
                final_pfx = dedup_prefix
        return (final_pfx, True)
    elif chr_filter:
        logging.info(f"Filtering '{stem_name}' to chromosome(s) {chr_filter}...")
        subset_prefix = f'{out_prefix}_chr'
        rc = run([plink_exe, '--bfile', converted, '--chr', chr_filter, '--snps-only', 'just-acgt', '--set-hh-missing', '--allow-no-sex', '--make-bed', '--out', subset_prefix], quiet_fail=True)
        if rc == 0 and Path(f'{subset_prefix}.bed').exists():
            if is_temp:
                for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                    Path(f'{converted}{ext}').unlink(missing_ok=True)
            return (subset_prefix, True)
    return (converted, is_temp)

def run_convert_pipeline(datasets: list[str], out_prefixes: list[str] | None, to_plink: bool, to_eig: bool, chr_filter: str | None, target_fids: set[str], target_iids: set[str], king_cutoff: float | None, default_fid: str | None, default_iid: str | None, sex: int, plink_exe: str, plink2_exe: str, geno: float | None=None, mind: float | None=None, maf: float | None=None, mac: int | None=None, exclude_high_ld: bool=False, ld_params: tuple[int, int, float] | None=None, extract_snps: str | None=None, exclude_snps: str | None=None) -> None:
    if to_plink and to_eig:
        sys.exit('ERROR: Cannot specify both --plink and --eig for conversion. Choose one target format.')
    total_ds = len(datasets)
    ui_step(1, 1, f'Running dataset conversion on {total_ds} dataset(s)...')
    for idx, src in enumerate(datasets, 1):
        stem = get_dataset_stem(src)
        fmt = detect_format(src)
        if fmt == 'error':
            p = Path(src)
            if p.name.lower().endswith(('.vcf', '.vcf.gz')):
                fmt = 'vcf'
            elif p.name.lower().endswith(('.bed', '.bim', '.fam')):
                fmt = 'plink'
            elif p.name.lower().endswith(('.geno', '.snp', '.ind')):
                fmt = 'eigenstrat'
            elif p.name.lower().endswith(('.txt', '.tsv', '.csv')):
                fmt = 'raw'
            else:
                sys.exit(f"ERROR: Cannot identify format for dataset '{src}'. Supported: PLINK binary, EIGENSTRAT, VCF, or raw DTC.")
        if to_plink:
            target_fmt = 'plink'
        elif to_eig:
            target_fmt = 'eigenstrat'
        elif fmt == 'eigenstrat':
            target_fmt = 'plink'
        elif fmt == 'plink':
            target_fmt = 'eigenstrat'
        else:
            target_fmt = 'plink'
        if out_prefixes and idx - 1 < len(out_prefixes):
            out_pfx = out_prefixes[idx - 1]
        elif out_prefixes and len(out_prefixes) == 1 and (total_ds == 1):
            out_pfx = out_prefixes[0]
        else:
            out_pfx = f'{stem}_plink' if target_fmt == 'plink' else f'{stem}_eig'
            if out_pfx == src:
                out_pfx = f'{stem}_converted'
        Path(out_pfx).parent.mkdir(parents=True, exist_ok=True)
        ui_info(f"[{idx}/{total_ds}] '{src}' ({fmt.upper()}) -> {target_fmt.upper()} (Output: '{out_pfx}')...")
        resolved_path = resolve_raw_path(src) or Path(src)
        temp_prefixes: list[str] = []
        is_temp = False
        if fmt == 'vcf':
            vcf_pfx = f'{out_pfx}_vcf_tmp'
            sfid = default_fid or stem
            siid = default_iid or stem
            plink_pfx = convert_vcf_to_plink(resolved_path, vcf_pfx, sfid, siid, sex, plink2_exe)
            is_temp = True
            temp_prefixes.append(plink_pfx)
        elif fmt in ('raw', '23andme'):
            raw_pfx = f'{out_pfx}_raw_tmp'
            sfid = default_fid or stem
            siid = default_iid or stem
            plink_pfx = convert_raw_dtc_to_plink(resolved_path, raw_pfx, sfid, siid, sex, plink_exe, plink2_exe)
            is_temp = True
            temp_prefixes.append(plink_pfx)
        elif fmt == 'eigenstrat':
            bk_col(src)
            nb_col(src)
            eig_pfx = f'{out_pfx}_eig2plink_tmp'
            plink_pfx = convert_to_plink(src, plink2_exe, out_prefix=eig_pfx)
            ps_col(plink_pfx, src)
            shutil.copy2(f'{src}_backup.ind', f'{src}.ind')
            Path(f'{src}_backup.ind').unlink(missing_ok=True)
            is_temp = True
            temp_prefixes.append(plink_pfx)
        elif fmt == 'plink':
            p = Path(src)
            plink_stem = str(p.parent / p.stem) if p.name.endswith(('.bed', '.bim', '.fam')) else src
            plink_pfx = plink_stem
            is_temp = False
        if target_fids or target_iids:
            filter_pfx = f'{out_pfx}_filtered'
            filtered_res, filter_is_temp = filter_plink_samples(bfile_prefix=plink_pfx, out_prefix=filter_pfx, plink_exe=plink_exe, target_fids=target_fids, target_iids=target_iids, dataset_name=f"Dataset '{src}'")
            if filter_is_temp:
                temp_prefixes.append(filtered_res)
                if is_temp:
                    for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                        Path(f'{plink_pfx}{ext}').unlink(missing_ok=True)
                plink_pfx = filtered_res
                is_temp = True
        if chr_filter:
            chr_pfx = f'{out_pfx}_chr'
            rc = run([plink_exe, '--bfile', plink_pfx, '--chr', chr_filter, '--allow-no-sex', '--make-bed', '--out', chr_pfx], quiet_fail=True)
            if rc == 0 and Path(f'{chr_pfx}.bed').exists():
                temp_prefixes.append(chr_pfx)
                if is_temp:
                    for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                        Path(f'{plink_pfx}{ext}').unlink(missing_ok=True)
                plink_pfx = chr_pfx
                is_temp = True
            else:
                cleanup(out_pfx, extra_prefixes=temp_prefixes)
                sys.exit(f"ERROR: Failed to filter chromosome(s) {chr_filter} for '{src}'")
        if king_cutoff is not None:
            king_pfx, removed_king, was_king = apply_king_cutoff(bfile_prefix=plink_pfx, out_prefix=f'{out_pfx}_king', cutoff=king_cutoff, plink2_exe=plink2_exe, plink_exe=plink_exe, dataset_label=f"Dataset '{src}'")
            if was_king:
                temp_prefixes.append(king_pfx)
                if is_temp:
                    for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                        Path(f'{plink_pfx}{ext}').unlink(missing_ok=True)
                plink_pfx = king_pfx
                is_temp = True
                ui_ok(f'KING filter ({king_cutoff}): removed {len(removed_king)} sample(s)')
        curr_qc, qc_temp = apply_qc_filters(bfile_prefix=plink_pfx, out_prefix=f'{out_pfx}_qc', plink_exe=plink_exe, geno=geno, mind=mind, maf=maf, mac=mac, exclude_high_ld=exclude_high_ld)
        if qc_temp:
            temp_prefixes.append(curr_qc)
            if is_temp:
                for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                    Path(f'{plink_pfx}{ext}').unlink(missing_ok=True)
            plink_pfx = curr_qc
            is_temp = True
        if ld_params is not None:
            ld_w, ld_s, ld_r = ld_params
            curr_ld, ld_temp = apply_ld_pruning(bfile_prefix=plink_pfx, out_prefix=f'{out_pfx}_ld', plink_exe=plink_exe, window=ld_w, step=ld_s, r2=ld_r)
            if ld_temp:
                temp_prefixes.append(curr_ld)
                if is_temp:
                    for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                        Path(f'{plink_pfx}{ext}').unlink(missing_ok=True)
                plink_pfx = curr_ld
                is_temp = True
        if extract_snps or exclude_snps:
            curr_vf, vf_temp = apply_variant_filters(bfile_prefix=plink_pfx, out_prefix=f'{out_pfx}_snps', plink_exe=plink_exe, extract_file=extract_snps, exclude_file=exclude_snps)
            if vf_temp:
                temp_prefixes.append(curr_vf)
                if is_temp:
                    for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                        Path(f'{plink_pfx}{ext}').unlink(missing_ok=True)
                plink_pfx = curr_vf
                is_temp = True
        if target_fmt == 'plink':
            if plink_pfx != out_pfx:
                for ext in ('.bed', '.bim', '.fam'):
                    s = Path(f'{plink_pfx}{ext}')
                    d = Path(f'{out_pfx}{ext}')
                    if is_temp:
                        shutil.move(s, d)
                    else:
                        shutil.copy2(s, d)
            format_plink_fam(Path(f'{out_pfx}.fam'))
            final_ns, final_nv = get_dataset_counts(out_pfx)
            cleanup(out_pfx, delete_final_plink=False, extra_prefixes=temp_prefixes)
            ui_ok(f"[{idx}/{total_ds}] Converted '{src}' -> PLINK binary ({final_ns:,} samples, {final_nv:,} variants):")
            print(f'       {out_pfx}.bed')
            print(f'       {out_pfx}.bim')
            print(f'       {out_pfx}.fam\n')
        elif target_fmt == 'eigenstrat':
            tmp_eig_pfx = f'{out_pfx}_eig_conv_tmp'
            for ext in ('.bed', '.bim', '.fam'):
                shutil.copy2(f'{plink_pfx}{ext}', f'{tmp_eig_pfx}{ext}')
            temp_prefixes.append(tmp_eig_pfx)
            if is_temp:
                for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                    Path(f'{plink_pfx}{ext}').unlink(missing_ok=True)
            export_eigenstrat(tmp_eig_pfx, plink2_exe)
            for ext in ('.geno', '.snp', '.ind'):
                shutil.move(Path(f'{tmp_eig_pfx}{ext}'), Path(f'{out_pfx}{ext}'))
            final_ns, final_nv = get_dataset_counts(out_pfx)
            cleanup(out_pfx, delete_final_plink=False, extra_prefixes=temp_prefixes)
            ui_ok(f"[{idx}/{total_ds}] Converted '{src}' -> EIGENSTRAT ({final_ns:,} samples, {final_nv:,} variants):")
            print(f'       {out_pfx}.geno')
            print(f'       {out_pfx}.snp')
            print(f'       {out_pfx}.ind\n')
    print(f'[DONE] Successfully converted {total_ds} dataset(s).\n')

def main() -> None:
    parser = argparse.ArgumentParser(description='Merge genomic datasets (PLINK binary, EIGENSTRAT, commercial raw arrays, or VCF) into a reference dataset.', formatter_class=argparse.RawDescriptionHelpFormatter, epilog='\nExamples:\n  # Standard single PLINK or EIGENSTRAT merge:\n  genomerge --ref Data_Ref --sub Data1 --out test_single\n\n  # Multi-cohort CLI merge:\n  genomerge --ref Data_Ref --sub cohort1 cohort2 --out test_multi\n\n  # Batch file merge:\n  genomerge --ref Data_Ref --batch list.txt --out test_batch\n\n  # Mixed formats:\n  genomerge --ref Data_Ref --sub Data1 --raw ancestry.txt --vcf wgs.vcf.gz --out test_mixed\n\n  # Autosomal-only merge:\n  genomerge --ref Data_Ref --sub Data1 --autosome --out test_autosome\n\n  # Autosomal merge with KING kinship duplicate removal:\n  genomerge --ref Data_Ref --sub Data1 --autosome --king-cutoff 0.354 --out test_king\n\n  # Standalone dataset conversion to EIGENSTRAT:\n  genomerge --convert Data1 --eig --out Data1_eig\n\n  # Standalone dataset conversion to PLINK binary:\n  genomerge --convert Data_Ref --plink --out Data_Ref_plink\n\n  # Multi-dataset batch conversion:\n  genomerge --convert Data1 Data2 Data3 --plink --out Data1_p Data2_p Data3_p\n\n  # Convert commercial DTC raw array or VCF with custom sample info:\n  genomerge --convert 23andme.txt --plink --out my_sample --pop Population1 --id Sample1\n\n  # Dry-run pre-merge integrity check:\n  genomerge --ref Data_Ref --sub Data1 --check\n        ')
    parser.add_argument('--check-tools', action='store_true', help='Check availability of external dependencies (PLINK 1.9 and PLINK 2.0) and exit')
    parser.add_argument('--ld-1', action='store_true', help='Modern LD pruning preset: --indep-pairwise 50 10 0.2')
    parser.add_argument('--ld-2', action='store_true', help='Ancient DNA (aDNA) LD pruning preset: --indep-pairwise 200 25 0.4')
    parser.add_argument('--ld', '--ld-prune', dest='ld_args', nargs='*', default=None, help="Custom LD pruning: specify window step r2 (e.g. --ld 50 10 0.2), or '1' (modern) / '2' (aDNA)")
    parser.add_argument('--exclude-high-ld', action='store_true', help='Exclude long-range high-LD inversion regions (HLA, 8p23, 17q21, LCT)')
    parser.add_argument('--maf', type=float, default=None, help='Filter out variants with MAF < threshold (e.g. 0.01)')
    parser.add_argument('--mac', type=int, default=None, help='Filter out variants with MAC < threshold (e.g. 2)')
    parser.add_argument('--geno', type=float, default=None, help='Filter out variants with missing rate > threshold (e.g. 0.10)')
    parser.add_argument('--mind', type=float, default=None, help='Filter out samples with missing rate > threshold (e.g. 0.20)')
    parser.add_argument('--extract', dest='extract_snps', default=None, metavar='FILE', help='Extract variants listed in a PLINK .snplist / text file (one variant ID per line)')
    parser.add_argument('--exclude', dest='exclude_snps', default=None, metavar='FILE', help='Exclude variants listed in a PLINK .snplist / text file (one variant ID per line)')
    parser.add_argument('--ref', required=False, default=None, help='Reference dataset prefix (PLINK or EIGENSTRAT)')
    parser.add_argument('--sub', nargs='*', default=None, help='Subject dataset prefix(es)')
    parser.add_argument('--convert', nargs='+', default=None, help='Convert one or more datasets to PLINK binary or EIGENSTRAT without merging')
    parser.add_argument('--raw', '--dtc', dest='raw_file', nargs='*', default=None, help='Raw commercial DNA array file(s) (AncestryDNA, MyHeritage, FTDNA, LivingDNA, 23andMe)')
    parser.add_argument('--23', dest='twenty_three', nargs='*', default=None, help='23andMe file(s) for subject dataset (alias for --raw)')
    parser.add_argument('--vcf', dest='vcf_file', nargs='*', default=None, help='VCF format file(s) (.vcf or .vcf.gz)')
    parser.add_argument('--batch', nargs='*', default=None, help='Batch file(s) listing subject datasets, one per line (format hint optional)')
    parser.add_argument('--out', nargs='*', default=None, help='Output dataset prefix(es) (for merge or conversion)')
    parser.add_argument('--check', action='store_true', help='Dry run pre-merge compatibility checks and exit without writing merged output')
    parser.add_argument('--eig', '--eigenstrat', action='store_true', help='Output dataset(s) as EIGENSTRAT (.geno/.snp/.ind)')
    parser.add_argument('--plink', '--bfile', action='store_true', help='Output dataset(s) as PLINK binary (.bed/.bim/.fam) instead of EIGENSTRAT')
    parser.add_argument('--autosome', '--autosome-only', action='store_true', help='Filter datasets to autosomal chromosomes (1-22) only')
    parser.add_argument('--chr', nargs='*', default=None, help='Filter to specific chromosome(s) (e.g. 1-22, 21,22 or 21 22)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Show full raw PLINK and subprocess command output in console (default: quiet)')
    parser.add_argument('--keep-pop', '--keep-ref-pop', '--poplist', '--ref-poplist', dest='keep_pop', nargs='*', default=None, help='Text file(s) or population name(s) (FID) to keep from reference')
    parser.add_argument('--keep-ind', '--keep-ref-ind', '--indlist', '--ref-indlist', dest='keep_ind', nargs='*', default=None, help='Text file(s) or sample ID(s) (IID) to keep from reference')
    parser.add_argument('--keep-sub-pop', '--sub-poplist', dest='keep_sub_pop', nargs='*', default=None, help='Text file(s) or population name(s) (FID) to keep from subject cohort(s)')
    parser.add_argument('--keep-sub-ind', '--sub-indlist', dest='keep_sub_ind', nargs='*', default=None, help='Text file(s) or sample ID(s) (IID) to keep from subject cohort(s)')
    parser.add_argument('--keep-all', action='store_true', help='Apply --keep-pop and --keep-ind to both reference and subject cohorts')
    parser.add_argument('--sub-pop', '--pop', default=None, help='Population label (FID) for subject sample (default: filename stem)')
    parser.add_argument('--sub-id', '--id', default=None, help='Sample ID (IID) for subject sample (default: filename stem)')
    parser.add_argument('--sub-sex', '--sex', type=int, default=0, help='Sex for subject sample: 1=male, 2=female, 0=unknown (default: 0)')
    parser.add_argument('--king-cutoff', '--king', dest='king_cutoff', type=float, nargs='?', const=0.354, default=None, help='Filter merged dataset using KING kinship cutoff (default: 0.354 for duplicate/twin removal; e.g. 0.177 for 1st-degree)')
    parser.add_argument('--king-ref', dest='king_ref', type=float, nargs='?', const=0.354, default=None, help='Filter reference dataset using KING kinship cutoff before merging (default: 0.354)')
    parser.add_argument('--king-sub', dest='king_sub', type=float, nargs='?', const=0.354, default=None, help='Filter subject dataset(s) using KING kinship cutoff before merging (default: 0.354)')
    parser.add_argument('--king-all', dest='king_all', type=float, nargs='?', const=0.354, default=None, help='Apply KING kinship cutoff to reference, subjects, and merged dataset (default: 0.354)')
    args = parser.parse_args()
    if args.extract_snps and not Path(args.extract_snps).is_file():
        sys.exit(f"ERROR: --extract file '{args.extract_snps}' not found")
    if args.exclude_snps and not Path(args.exclude_snps).is_file():
        sys.exit(f"ERROR: --exclude file '{args.exclude_snps}' not found")
    script_dir = Path(__file__).resolve().parent
    if args.check_tools:
        ok = check_tools(script_dir)
        sys.exit(0 if ok else 1)
    ld_params = None
    if args.ld_1:
        ld_params = (50, 10, 0.2)
    elif args.ld_2:
        ld_params = (200, 25, 0.4)
    elif args.ld_args is not None:
        if len(args.ld_args) == 0 or (len(args.ld_args) == 1 and args.ld_args[0] in ('1', 'modern')):
            ld_params = (50, 10, 0.2)
        elif len(args.ld_args) == 1 and args.ld_args[0] in ('2', 'adna', 'ancient'):
            ld_params = (200, 25, 0.4)
        elif len(args.ld_args) >= 3:
            try:
                ld_params = (int(args.ld_args[0]), int(args.ld_args[1]), float(args.ld_args[2]))
            except ValueError:
                sys.exit(f"ERROR: Invalid arguments for --ld: expected '<window:int> <step:int> <r2:float>', got {args.ld_args}")
        else:
            sys.exit('ERROR: --ld requires 3 values: <window> <step> <r2>, or use --ld-1 / --ld-2')
    if args.chr:
        chr_filter = ','.join((c.strip() for c in args.chr if c.strip()))
    elif args.autosome:
        chr_filter = '1-22'
    else:
        chr_filter = None
    default_king = args.king_all if args.king_all is not None else args.king_cutoff
    king_cutoff = args.king_cutoff if args.king_cutoff is not None else default_king
    king_ref = args.king_ref if args.king_ref is not None else default_king
    king_sub = args.king_sub if args.king_sub is not None else default_king
    script_dir = Path(__file__).resolve().parent
    plink_exe = find_exe('plink', script_dir)
    plink2_exe = find_exe('plink2', script_dir)
    ref_target_fids = load_id_list(args.keep_pop)
    ref_target_iids = load_id_list(args.keep_ind)
    if args.convert:
        log_out = args.out[0] if args.out and len(args.out) > 0 else 'convert'
        setup_logging(log_out, to_file=True, verbose=args.verbose)
        logging.info(f'plink:  {plink_exe}')
        logging.info(f'plink2: {plink2_exe}')
        run_convert_pipeline(datasets=args.convert, out_prefixes=args.out, to_plink=args.plink, to_eig=args.eig, chr_filter=chr_filter, target_fids=ref_target_fids | load_id_list(args.keep_sub_pop), target_iids=ref_target_iids | load_id_list(args.keep_sub_ind), king_cutoff=king_cutoff, default_fid=args.sub_pop, default_iid=args.sub_id, sex=args.sub_sex, plink_exe=plink_exe, plink2_exe=plink2_exe, geno=args.geno, mind=args.mind, maf=args.maf, mac=args.mac, exclude_high_ld=args.exclude_high_ld, ld_params=ld_params, extract_snps=args.extract_snps, exclude_snps=args.exclude_snps)
        return
    if not args.ref:
        parser.error('the following arguments are required: --ref (or use --convert for standalone dataset conversion)')
    if not args.check and (not args.out):
        parser.error('the following arguments are required: --out')
    out_prefix = args.out[0] if isinstance(args.out, list) and args.out else args.out or 'check_temp'
    args.out = out_prefix
    write_log_file = args.out is not None or not args.check
    setup_logging(out_prefix, to_file=write_log_file, verbose=args.verbose)
    logging.info(f'plink:  {plink_exe}')
    logging.info(f'plink2: {plink2_exe}')
    if args.keep_all:
        sub_target_fids = set(ref_target_fids)
        sub_target_iids = set(ref_target_iids)
    else:
        sub_target_fids = load_id_list(args.keep_sub_pop)
        sub_target_iids = load_id_list(args.keep_sub_ind)
    if ref_target_fids:
        logging.info(f'Loaded {len(ref_target_fids)} target populations (FID) to keep in reference')
    if ref_target_iids:
        logging.info(f'Loaded {len(ref_target_iids)} target individuals (IID) to keep in reference')
    if sub_target_fids:
        logging.info(f'Loaded {len(sub_target_fids)} target populations (FID) to keep in subjects')
    if sub_target_iids:
        logging.info(f'Loaded {len(sub_target_iids)} target individuals (IID) to keep in subjects')
    raw_flag = args.raw_file == [] or args.twenty_three == []
    vcf_flag = args.vcf_file == []
    raw_subject_entries: list[tuple[str, str | None]] = []
    if args.sub is not None:
        for s in args.sub:
            if raw_flag:
                raw_subject_entries.append((s, 'raw'))
            elif vcf_flag:
                raw_subject_entries.append((s, 'vcf'))
            else:
                raw_subject_entries.append((s, None))
    if args.raw_file and args.raw_file != []:
        for r in args.raw_file:
            raw_subject_entries.append((r, 'raw'))
    if args.twenty_three and args.twenty_three != []:
        for r in args.twenty_three:
            raw_subject_entries.append((r, 'raw'))
    if args.vcf_file and args.vcf_file != []:
        for v in args.vcf_file:
            raw_subject_entries.append((v, 'vcf'))
    if args.batch is not None:
        for b in args.batch:
            parsed = parse_batch_file(Path(b))
            raw_subject_entries.extend(parsed)
    if not raw_subject_entries:
        sys.exit('ERROR: Must specify at least one subject dataset using --sub, --batch, --raw, or --vcf')
    temp_prefixes: list[str] = []
    total_phases = 3 if args.check else 4 if args.plink else 5
    ui_step(1, total_phases, f"Loading Reference: '{args.ref}'")
    ref_plink, ref_is_temp = normalize_ref_to_plink(args.ref, out_prefix, plink2_exe)
    if ref_is_temp:
        temp_prefixes.append(ref_plink)
    ref_ns, ref_nv = get_dataset_counts(ref_plink)
    if ref_target_fids or ref_target_iids:
        filter_pfx = f'{out_prefix}_ref_filtered'
        ref_filtered, filtered_is_temp = filter_plink_samples(bfile_prefix=ref_plink, out_prefix=filter_pfx, plink_exe=plink_exe, target_fids=ref_target_fids, target_iids=ref_target_iids, dataset_name=f"Reference '{args.ref}'")
        if filtered_is_temp:
            temp_prefixes.append(ref_filtered)
            ref_plink = ref_filtered
            ref_ns, ref_nv = get_dataset_counts(ref_plink)
    if chr_filter:
        need_ref_filter = True
        if chr_filter == '1-22' and (not bim_has_other_chromosomes(ref_plink, AUTOSOMAL_CHRS)):
            need_ref_filter = False
        if need_ref_filter:
            logging.info(f"Filtering reference '{ref_plink}' to chromosome(s) {chr_filter}...")
            ref_chr_pfx = f'{out_prefix}_ref_chr'
            rc = run([plink_exe, '--bfile', ref_plink, '--chr', chr_filter, '--allow-no-sex', '--make-bed', '--out', ref_chr_pfx], quiet_fail=True)
            if rc == 0 and Path(f'{ref_chr_pfx}.bed').exists():
                temp_prefixes.append(ref_chr_pfx)
                ref_plink = ref_chr_pfx
                ref_ns, ref_nv = get_dataset_counts(ref_plink)
            else:
                cleanup(out_prefix, extra_prefixes=temp_prefixes)
                sys.exit(f"ERROR: Failed to filter reference '{ref_plink}' to chromosome(s) {chr_filter}")
    if king_ref is not None:
        ui_info(f"Filtering Reference '{args.ref}' using KING kinship cutoff (threshold: {king_ref})...")
        ref_king, removed_ref, was_ref_filtered = apply_king_cutoff(bfile_prefix=ref_plink, out_prefix=f'{out_prefix}_ref_king', cutoff=king_ref, plink2_exe=plink2_exe, plink_exe=plink_exe, dataset_label=f"Reference '{args.ref}'")
        if was_ref_filtered:
            temp_prefixes.append(ref_king)
            ref_plink = ref_king
            ref_ns, ref_nv = get_dataset_counts(ref_plink)
            ui_ok(f'Reference KING filter ({king_ref}): removed {len(removed_ref)} related sample(s), kept {ref_ns:,} samples')
            logging.info(f'Reference KING filter ({king_ref}) removed {len(removed_ref)} sample(s): {removed_ref}')
        else:
            ui_ok(f'Reference KING filter ({king_ref}): zero related pairs detected (> {king_ref})')
    ref_snplist_pfx = f'{out_prefix}_ref_snplist'
    snplist_cmd = [plink_exe, '--bfile', ref_plink]
    if chr_filter:
        snplist_cmd.extend(['--chr', chr_filter])
    snplist_cmd.extend(['--snps-only', 'just-acgt', '--write-snplist', '--allow-no-sex', '--out', ref_snplist_pfx])
    rc = run(snplist_cmd)
    ref_snplist_file = Path(f'{ref_snplist_pfx}.snplist')
    if rc != 0 or not ref_snplist_file.is_file():
        cleanup(out_prefix, extra_prefixes=temp_prefixes)
        sys.exit(f"ERROR: Failed to extract reference SNP list from '{ref_plink}'")
    temp_prefixes.append(ref_snplist_pfx)
    ui_ok(f'Reference ready ({ref_ns:,} samples, {ref_nv:,} variants)')
    ui_step(2, total_phases, f'Normalizing {len(raw_subject_entries)} Subject Cohort(s) to Ref SNP subset...')
    normalized_subjects: list[tuple[str, str, bool]] = []
    is_multi = len(raw_subject_entries) > 1
    for idx, (sub_raw, fmt_hint) in enumerate(raw_subject_entries, 1):
        stem = get_dataset_stem(sub_raw)
        if is_multi:
            default_fid = args.sub_pop or stem
            default_iid = args.sub_id or stem
        else:
            default_fid = args.sub_pop
            default_iid = args.sub_id
        norm_prefix = f'{out_prefix}_Sub_{idx}_norm'
        sub_plink, sub_is_temp = normalize_subject_to_plink(sub_raw=sub_raw, format_hint=fmt_hint, out_prefix=norm_prefix, default_fid=default_fid, default_iid=default_iid, sex=args.sub_sex, plink_exe=plink_exe, plink2_exe=plink2_exe, ref_snplist_path=ref_snplist_file, chr_filter=chr_filter)
        if sub_target_fids or sub_target_iids:
            filter_sub_pfx = f'{norm_prefix}_filtered'
            sub_filtered, filter_sub_is_temp = filter_plink_samples(bfile_prefix=sub_plink, out_prefix=filter_sub_pfx, plink_exe=plink_exe, target_fids=sub_target_fids, target_iids=sub_target_iids, dataset_name=f"Subject '{stem}'", allow_empty=True)
            if sub_filtered is None:
                ui_warn(f"Subject '{stem}' has 0 samples matching subject keep criteria. Skipping.")
                if sub_is_temp:
                    for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                        Path(f'{sub_plink}{ext}').unlink(missing_ok=True)
                continue
            if filter_sub_is_temp:
                if sub_is_temp:
                    for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                        Path(f'{sub_plink}{ext}').unlink(missing_ok=True)
                sub_plink = sub_filtered
                sub_is_temp = True
        if king_sub is not None:
            sub_king, removed_sub, was_sub_filtered = apply_king_cutoff(bfile_prefix=sub_plink, out_prefix=f'{norm_prefix}_king', cutoff=king_sub, plink2_exe=plink2_exe, plink_exe=plink_exe, dataset_label=f"Subject '{stem}'")
            if was_sub_filtered:
                if sub_is_temp:
                    for ext in ('.bed', '.bim', '.fam', '.log', '.nosex'):
                        Path(f'{sub_plink}{ext}').unlink(missing_ok=True)
                sub_plink = sub_king
                sub_is_temp = True
                temp_prefixes.append(sub_plink)
                sub_ns, sub_nv = get_dataset_counts(sub_plink)
                ui_ok(f"[{idx}/{len(raw_subject_entries)}] '{stem}': KING filter ({king_sub}) removed {len(removed_sub)} sample(s), kept {sub_ns:,} samples")
            else:
                ui_ok(f"[{idx}/{len(raw_subject_entries)}] '{stem}': KING filter ({king_sub}) zero related pairs (> {king_sub})")
        if sub_is_temp:
            temp_prefixes.append(sub_plink)
        normalized_subjects.append((sub_raw, sub_plink, sub_is_temp))
        sub_ns, sub_nv = get_dataset_counts(sub_plink)
        ui_ok(f"[{idx}/{len(raw_subject_entries)}] '{stem}': ready ({sub_ns:,} samples, {sub_nv:,} variants [subset of ref])")
    ui_step(3, total_phases, 'Pre-Merge Validation & Compatibility Checks...')
    all_datasets_for_check = [(f'Ref({args.ref})', ref_plink)] + [(f'Sub({lbl})', pfx) for lbl, pfx, _ in normalized_subjects]
    duplicates = check_duplicate_samples(all_datasets_for_check)
    if duplicates:
        cleanup(out_prefix, extra_prefixes=temp_prefixes)
        err_lines = ['Duplicate sample IDs (FID, IID) detected across datasets:']
        for (fid, iid), labels in duplicates:
            err_lines.append(f"  - Sample FID='{fid}', IID='{iid}' found in: {', '.join(labels)}")
        err_msg = '\n'.join(err_lines)
        ui_err(err_msg)
        logging.error(f'ERROR: {err_msg}')
        sys.exit(f'ERROR: {err_msg}')
    ui_ok('Sample ID collision check: OK (zero duplicates across cohorts)')
    ref_variants = get_bim_variants(ref_plink)
    valid_subjects: list[tuple[str, str, bool]] = []
    for lbl, pfx, is_tmp in normalized_subjects:
        sub_variants = get_bim_variants(pfx)
        overlap = len(ref_variants & sub_variants)
        if overlap == 0:
            ui_warn(f"Subject '{lbl}' has 0 variant overlap with Reference '{args.ref}'. Skipping.")
            logging.warning(f"WARNING: Subject '{lbl}' has 0 variant overlap with Reference '{args.ref}'. Skipping.")
            if is_tmp:
                cleanup_step(pfx, pfx, temp_sub=pfx)
        else:
            logging.info(f"Subject '{lbl}' has {overlap} overlapping variants with Reference '{args.ref}'.")
            valid_subjects.append((lbl, pfx, is_tmp))
    ui_ok(f'Variant overlap check: OK ({len(valid_subjects)}/{len(normalized_subjects)} cohorts valid)')
    if king_cutoff is not None or king_sub is not None or king_ref is not None:
        for lbl, pfx, _ in valid_subjects:
            c = king_sub or king_cutoff
            if c is not None:
                stem_lbl = get_dataset_stem(lbl)
                removed_sub_chk = check_king_cutoff(pfx, c, plink2_exe, f'{out_prefix}_chk_{stem_lbl}')
                if removed_sub_chk:
                    sample_summary = ', '.join((iid for _, iid in removed_sub_chk[:5]))
                    if len(removed_sub_chk) > 5:
                        sample_summary += f' ... (+{len(removed_sub_chk) - 5} more)'
                    ui_warn(f"KING check ({c}) on Subject '{stem_lbl}': {len(removed_sub_chk)} sample(s) flagged for removal: {sample_summary}")
                    logging.warning(f"KING check ({c}) on Subject '{stem_lbl}': flagged {removed_sub_chk}")
                else:
                    ui_ok(f"KING check ({c}) on Subject '{stem_lbl}': OK (zero related pairs > {c})")
        if king_ref is not None:
            removed_ref_chk = check_king_cutoff(ref_plink, king_ref, plink2_exe, f'{out_prefix}_chk_ref')
            if removed_ref_chk:
                sample_summary = ', '.join((iid for _, iid in removed_ref_chk[:5]))
                if len(removed_ref_chk) > 5:
                    sample_summary += f' ... (+{len(removed_ref_chk) - 5} more)'
                ui_warn(f'KING check ({king_ref}) on Reference: {len(removed_ref_chk)} sample(s) flagged for removal: {sample_summary}')
                logging.warning(f'KING check ({king_ref}) on Reference: flagged {removed_ref_chk}')
            else:
                ui_ok(f'KING check ({king_ref}) on Reference: OK (zero related pairs > {king_ref})')
    if args.check:
        print(f'\n[CHECK COMPLETED] Pre-merge checks passed successfully for {len(valid_subjects)} cohort(s).')
        logging.info('Pre-merge validation checks completed successfully. (--check dry run)')
        cleanup(out_prefix, extra_prefixes=temp_prefixes)
        sys.exit(0)
    if not valid_subjects:
        cleanup(out_prefix, extra_prefixes=temp_prefixes)
        ui_err('No valid subject datasets remaining after variant overlap check.')
        sys.exit('ERROR: No valid subject datasets remaining after variant overlap check.')
    num_subs = len(valid_subjects)
    curr_ref = ref_plink
    ui_step(4, total_phases, f'Merging {num_subs} Cohort(s) into Reference...')
    for idx, (sub_lbl, sub_pfx, is_tmp) in enumerate(valid_subjects, 1):
        is_last = idx == num_subs
        step_out = args.out if is_last else f'{args.out}_step{idx}'
        step_prefix = args.out if num_subs == 1 else f'{args.out}_step{idx}'
        ui_info(f"[{idx}/{num_subs}] Merging '{sub_lbl}' into reference...")
        logging.info(f"--- Cohort Merge Step {idx}/{num_subs}: Merging '{sub_lbl}' into reference ---")
        ok = merge_pipeline(ref=curr_ref, sub=sub_pfx, out=step_out, plink_exe=plink_exe, plink2_exe=plink2_exe, step_prefix=step_prefix, chr_filter=chr_filter)
        if not ok:
            cleanup(args.out, extra_prefixes=temp_prefixes)
            sys.exit(f"ERROR: Merge failed for subject '{sub_lbl}' (Step {idx})")
        prev_step_to_clean = curr_ref if idx > 1 else None
        temp_sub_to_clean = sub_pfx if is_tmp else None
        cleanup_step(step_prefix=step_prefix, step_out=step_out, prev_ref=prev_step_to_clean, temp_sub=temp_sub_to_clean)
        curr_ref = step_out
    if ref_is_temp:
        cleanup_step(ref_plink, ref_plink, prev_ref=ref_plink)
    if king_cutoff is not None:
        ui_info(f'Filtering merged cohort using KING kinship cutoff (threshold: {king_cutoff})...')
        filtered_king, removed_king, was_king_filtered = apply_king_cutoff(bfile_prefix=curr_ref, out_prefix=f'{args.out}_merged_king', cutoff=king_cutoff, plink2_exe=plink2_exe, plink_exe=plink_exe, dataset_label='merged cohort')
        if was_king_filtered:
            temp_prefixes.append(filtered_king)
            for ext in ('.bed', '.bim', '.fam'):
                src = Path(f'{filtered_king}{ext}')
                dst = Path(f'{args.out}{ext}')
                if src.exists():
                    shutil.copy2(src, dst)
            final_kns, final_knv = get_dataset_counts(args.out)
            sample_summary = ', '.join((iid for _, iid in removed_king[:5]))
            if len(removed_king) > 5:
                sample_summary += f' ... (+{len(removed_king) - 5} more)'
            ui_ok(f'Merged KING filter ({king_cutoff}): removed {len(removed_king)} sample(s) ({sample_summary}), kept {final_kns:,} samples')
            logging.info(f'Merged KING filter ({king_cutoff}): removed {len(removed_king)} sample(s): {removed_king}')
        else:
            ui_ok(f'Merged KING filter ({king_cutoff}): OK (zero related pairs > {king_cutoff})')
    curr_ref, qc_temp = apply_qc_filters(bfile_prefix=curr_ref, out_prefix=f'{args.out}_merged_qc', plink_exe=plink_exe, geno=args.geno, mind=args.mind, maf=args.maf, mac=args.mac, exclude_high_ld=args.exclude_high_ld)
    if qc_temp:
        temp_prefixes.append(curr_ref)
        for ext in ('.bed', '.bim', '.fam'):
            shutil.copy2(f'{curr_ref}{ext}', f'{args.out}{ext}')
        curr_ref = args.out
    if ld_params is not None:
        ld_w, ld_s, ld_r = ld_params
        curr_ref, ld_temp = apply_ld_pruning(bfile_prefix=curr_ref, out_prefix=f'{args.out}_merged_ld', plink_exe=plink_exe, window=ld_w, step=ld_s, r2=ld_r)
        if ld_temp:
            temp_prefixes.append(curr_ref)
            for ext in ('.bed', '.bim', '.fam'):
                shutil.copy2(f'{curr_ref}{ext}', f'{args.out}{ext}')
            curr_ref = args.out
    if args.extract_snps or args.exclude_snps:
        curr_ref, vf_temp = apply_variant_filters(bfile_prefix=curr_ref, out_prefix=f'{args.out}_merged_snps', plink_exe=plink_exe, extract_file=args.extract_snps, exclude_file=args.exclude_snps)
        if vf_temp:
            temp_prefixes.append(curr_ref)
            for ext in ('.bed', '.bim', '.fam'):
                shutil.copy2(f'{curr_ref}{ext}', f'{args.out}{ext}')
            curr_ref = args.out
    if args.plink:
        format_plink_fam(Path(f'{args.out}.fam'))
        final_ns, final_nv = get_dataset_counts(args.out)
        cleanup(args.out, delete_final_plink=False, extra_prefixes=temp_prefixes)
        ui_ok(f'Validated and formatted {args.out}.fam (6 standard columns)')
        print(f'\n[DONE] Merged PLINK binary complete ({final_ns:,} samples, {final_nv:,} variants):')
        print(f'       {args.out}.bed')
        print(f'       {args.out}.bim')
        print(f'       {args.out}.fam')
        if write_log_file:
            print(f'       Full log: {out_prefix}_merge.log')
        logging.info(f'DONE. Output: {args.out}.bed / .bim / .fam')
    else:
        ui_step(5, total_phases, 'Exporting merged dataset to EIGENSTRAT format...')
        export_eigenstrat(args.out, plink2_exe)
        final_ns, final_nv = get_dataset_counts(args.out)
        cleanup(args.out, delete_final_plink=True, extra_prefixes=temp_prefixes)
        ui_ok(f'Exported {args.out}.geno / .snp / .ind')
        print(f'\n[DONE] Merged dataset complete ({final_ns:,} samples, {final_nv:,} variants):')
        print(f'       {args.out}.geno')
        print(f'       {args.out}.snp')
        print(f'       {args.out}.ind')
        if write_log_file:
            print(f'       Full log: {out_prefix}_merge.log')
        logging.info(f'DONE. Output: {args.out}.geno / .snp / .ind')
if __name__ == '__main__':
    main()
