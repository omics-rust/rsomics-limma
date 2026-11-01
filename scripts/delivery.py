"""Package, install and archive the scoped native limma CLI."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import tarfile

PRODUCT = 'rsomics-limma'
ROOT = Path(__file__).resolve().parents[1]


def identity(path):
    body = Path(path).read_bytes()
    return dict(bytes=len(body), sha256=hashlib.sha256(body).hexdigest())


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def package_contract(package):
    parent = Path(package['manifest_path']).parent
    result = {k:v for k,v in package.items() if k not in ('id', 'manifest_path', 'targets')}
    result['targets'] = [{**t, 'src_path': str(Path(t['src_path']).relative_to(parent))} for t in package['targets']]
    return result


def extract_package(source, destination, expected, head):
    destination = Path(destination)
    destination.mkdir(parents=True)
    with tarfile.open(source, 'r:gz') as archive:
        members = archive.getmembers()
        require(len(members) == len({m.name for m in members}), 'duplicate source package member')
        bodies = {}
        for member in members:
            path = PurePosixPath(member.name)
            require(not path.is_absolute() and '..' not in path.parts and path.parts[0] == destination.name,
                'unsafe source package path')
            require(member.isfile() or member.isdir(), 'unsupported source package member')
            if member.isfile():
                name = '/'.join(path.parts[1:])
                require(name and member.mode == 0o644, 'source package file mode differs')
                bodies[name] = archive.extractfile(member).read()
        require(set(bodies) == set(expected) | {'Cargo.toml.orig', '.cargo_vcs_info.json'},
            'complete source package membership differs')
        require(bodies['Cargo.toml.orig'] == expected['Cargo.toml'], 'original Cargo manifest differs')
        vcs = json.loads(bodies['.cargo_vcs_info.json'])
        require(vcs['git']['sha1'] == head and not vcs['git'].get('dirty', False) and vcs['path_in_vcs'] == '',
            'source package is not the exact clean CI head')
        for name, body in expected.items():
            if name != 'Cargo.toml':
                require(bodies[name] == body, 'source package body differs: ' + name)
        for name, body in bodies.items():
            path = destination / name
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open('xb') as stream:
                stream.write(body)
            path.chmod(0o644)
    return dict(**identity(source), members=len(bodies), original_members=len(expected),
        original_bodies_modes_match=True, source_commit=head, clean_vcs=True)


def source():
    temporary = Path(os.environ['RUNNER_TEMP'])
    delivery = temporary / 'delivery'
    metadata = json.loads((delivery / 'package-metadata.json').read_bytes())
    package, = [p for p in metadata['packages'] if p['name'] == PRODUCT]
    prefix = PRODUCT + '-' + package['version']
    source = Path(metadata['target_directory']) / 'package' / (prefix + '.crate')
    shutil.copy2(source, delivery / source.name)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    require(head == os.environ['GITHUB_SHA'], 'source head differs from actual CI head')
    require(subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT) == b'', 'source worktree is dirty')
    expected = {}
    for row in subprocess.check_output(['git', 'ls-tree', '-rz', head], cwd=ROOT).rstrip(b'\0').split(b'\0'):
        header, name = row.split(b'\t')
        mode, kind, blob = header.decode().split()
        name = name.decode()
        require(kind == 'blob' and mode == '100644', 'unexpected source Git type/mode: ' + name)
        body = (ROOT / name).read_bytes()
        require(hashlib.sha1(b'blob ' + str(len(body)).encode() + b'\0' + body).hexdigest() == blob,
            'source differs from exact Git body: ' + name)
        expected[name] = body
    destination = temporary / 'packaged-source' / prefix
    verification = extract_package(source, destination, expected, head)
    unpacked = json.loads(subprocess.check_output(['cargo', 'metadata', '--manifest-path', str(destination / 'Cargo.toml'),
        '--locked', '--offline', '--no-deps', '--format-version', '1']))
    installed, = [p for p in unpacked['packages'] if p['name'] == PRODUCT]
    require(package_contract(installed) == package_contract(package), 'normalized source package metadata differs')
    verification['normalized_metadata_match'] = True
    save(delivery / 'package-verification.json', verification)
    metadata_body = subprocess.check_output(['cargo', 'metadata', '--locked', '--format-version', '1'], cwd=ROOT)
    tree_body = subprocess.check_output(['cargo', 'tree', '--locked', '--target', os.environ['EXPECTED_HOST'],
        '--edges', 'normal,build', '--prefix', 'none', '--format', '{p}'], cwd=ROOT)
    (delivery / 'native-dependencies.metadata.json').write_bytes(metadata_body)
    (delivery / 'native-dependencies.tree.txt').write_bytes(tree_body)
    by_name = {(p['name'], p['version']): p for p in json.loads(metadata_body)['packages']}
    enabled = set()
    for line in tree_body.decode().splitlines():
        name, version = re.match(r'([^ ]+) v([^ ]+)', line).groups()
        if name != PRODUCT:
            enabled.add(by_name[name, version]['id'])
    inventory = json.loads((ROOT / 'LICENSES/provenance.json').read_bytes())
    require(enabled == set(inventory['target_graphs'][os.environ['EXPECTED_HOST']]), 'enabled target dependency inventory differs')
    require(identity(ROOT / 'Cargo.lock')['sha256'] == inventory['lock_sha256'], 'notice inventory lock differs')
    for row in inventory['copied_texts'].values():
        require(identity(ROOT / row['destination']) == {k:row[k] for k in ('bytes', 'sha256')}, 'bundled notice body differs')
    save(delivery / 'native-dependencies.json', dict(target=os.environ['EXPECTED_HOST'], enabled_packages=sorted(enabled),
        inventory_sha256=identity(ROOT / 'LICENSES/provenance.json')['sha256'],
        metadata_sha256=hashlib.sha256(metadata_body).hexdigest(), tree_sha256=hashlib.sha256(tree_body).hexdigest(),
        enabled_normal_build_graph_verified=True, linked_binary_inventory=False))
    with Path(os.environ['GITHUB_ENV']).open('a') as stream:
        stream.write('PACKAGE_DIRECTORY=' + str(destination) + '\n')
        stream.write('PACKAGE_VERSION=' + package['version'] + '\n')
        stream.write('SOURCE_PACKAGE=' + str(delivery / source.name) + '\n')


def table(path):
    with Path(path).open() as source:
        return [row for row in csv.reader(source, delimiter='\t') if row]


def check_lifecycle(directory, fixture, mode, bundle_identity):
    directory, fixture = Path(directory), Path(fixture)
    bundle = directory / (mode + '.json')
    require(identity(bundle) == bundle_identity, 'saved bundle changed during reload')
    analysis = json.loads(bundle.read_bytes())
    expected = table(fixture / ('counts.tsv' if mode == 'counts' else 'expression.tsv'))
    field = 'counts' if mode == 'counts' else 'expression'
    data = analysis['dataset'][field]
    genes, samples = [r[0] for r in expected[1:]], expected[0][1:]
    require(analysis['product'] == PRODUCT and analysis['workers'] == 4, 'installed product or workers differ')
    require(data['gene_ids'] == genes and data['sample_ids'] == samples and data['values'] == [float(v) for r in expected[1:] for v in r[1:]],
        'installed input identities or values differ')
    sample_rows = table(fixture / 'samples.tsv')
    sample_map = {r[0]:r[1:] for r in sample_rows[1:] if r}
    require(analysis['dataset']['metadata_columns'] == sample_rows[0][1:]
        and analysis['dataset']['metadata'] == [sample_map[sample] for sample in samples], 'installed sample metadata join differs')
    if mode == 'counts':
        require(analysis['schema'] == 1 and len(genes) == 13 and analysis['normalization']['kept_indices'] == list(range(11)),
            'installed count/filter contract differs')
        genes = genes[:11]
    else:
        require(analysis['schema'] == 2 and analysis['input_kind'] == 'log-expression' and len(genes) == 6
            and 'normalization' not in analysis, 'installed ordinary contract differs')
    state = json.loads((directory / (mode + '-results.json')).read_bytes())
    require(state == analysis['initial_results'] and state['full_base_reused'] and not state['linear_model_refit']
        and len(state['rows']) == len(genes), 'installed full saved-result reuse differs')
    results = table(directory / (mode + '-results.tsv'))
    require(results[0] == ['gene', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']
        and [r[0] for r in results[1:]] == genes and all(len(r) == 7 for r in results[1:]),
        'installed result table identities or shape differ')
    fields = ['log_fc', 'ave_expr', 't', 'p_value', 'fdr', 'lods']
    number = lambda value: value if value in ('NA', 'NaN', 'Inf', '-Inf') else float(value)
    require(all([number(v) for v in row[1:]] == [number(result[key]) for key in fields]
        for row, result in zip(results[1:], state['rows'])), 'installed result table differs from saved state')
    for kind in (('expression', 'weights', 'fitted', 'residuals') if mode == 'counts' else ('expression', 'fitted', 'residuals')):
        rows = table(directory / (mode + '-' + kind + '.tsv'))
        require(rows[0] == ['gene', *samples] and [r[0] for r in rows[1:]] == genes
            and all(len(r) == len(samples) + 1 for r in rows[1:]), 'installed export identities or shape differ: ' + kind)
        values = [float(v) for r in rows[1:] for v in r[1:]]
        require(all(float('-inf') < value < float('inf') for value in values), 'installed fixture export is not finite')
        if mode == 'expression' and kind == 'expression':
            require(values == data['values'], 'installed supplied-expression export changed')
    return dict(schema=analysis['schema'], input_features=len(expected)-1, output_features=len(genes), samples=len(samples),
        workers=4, full_saved_result_reused=True, bundle_unchanged=True)


def smoke(binary):
    binary = Path(binary).resolve(strict=True)
    directory = Path(os.environ['RUNNER_TEMP']) / 'delivery/lifecycle'
    directory.mkdir(parents=True)
    shutil.copy2(binary, directory / 'installed-binary')
    save(directory / 'installed-binary.json', identity(binary))
    fixture = ROOT / 'tests/data'
    records = []

    def run(arguments):
        command = [str(binary), *map(str, arguments)]
        row = dict(argv=command)
        records.append(row)
        save(directory / 'commands.json', records)
        try:
            with (directory / f'{len(records)-1}.stdout').open('xb') as out, (directory / f'{len(records)-1}.stderr').open('xb') as err:
                result = subprocess.run(command, cwd=directory, stdout=out, stderr=err, timeout=120)
            row['returncode'] = result.returncode
        except BaseException as error:
            row['error'] = repr(error)
            raise
        finally:
            save(directory / 'commands.json', records)
        require(result.returncode == 0, 'installed command failed: ' + str(command))

    run(['--version'])
    require((directory / '0.stdout').read_text().strip() == PRODUCT + ' ' + os.environ['PACKAGE_VERSION'], 'installed version differs')
    run(['--help'])
    contrast = ['--contrast=(Intercept)=0', '--contrast=typesingle-read=0', '--contrast=conditiontreated=1']
    checked = {}
    for mode in ('counts', 'expression'):
        bundle = directory / (mode + '.json')
        run(['--threads', '4', 'run', '--' + mode, fixture / (mode + '.tsv'), '--samples', fixture / 'samples.tsv',
            '--design', fixture / 'design.tsv', '--bundle', bundle, *contrast])
        before = identity(bundle)
        run(['results', '--bundle', bundle, '--output', directory / (mode + '-results.tsv'),
            '--state', directory / (mode + '-results.json'), *contrast])
        for kind in (('expression', 'weights', 'fitted', 'residuals') if mode == 'counts' else ('expression', 'fitted', 'residuals')):
            run(['export', '--bundle', bundle, '--kind', kind, '--output', directory / (mode + '-' + kind + '.tsv')])
        checked[mode] = check_lifecycle(directory, fixture, mode, before)
    save(directory / 'VERIFIED.json', dict(installed_binary=identity(binary), calls=records, workflows=checked,
        files={p.name:identity(p) for p in directory.iterdir() if p.is_file()},
        scope='installed package lifecycle only; not a new source or performance campaign'))


def bundle(binary):
    temporary = Path(os.environ['RUNNER_TEMP'])
    delivery = temporary / 'delivery'
    binary = Path(binary).resolve(strict=True)
    stem = PRODUCT + '-' + os.environ['PACKAGE_VERSION'] + '-' + os.environ['NATIVE_TARGET']
    folder = delivery / stem
    folder.mkdir()
    shutil.copy2(binary, folder / PRODUCT)
    for name in ('README.md', 'LICENSE-MIT', 'LICENSE-APACHE-2.0', 'THIRD_PARTY_NOTICES.md'):
        shutil.copy2(ROOT / name, folder / name)
    shutil.copytree(ROOT / 'LICENSES', folder / 'LICENSES')
    (folder / 'docs').mkdir()
    shutil.copy2(ROOT / 'docs/merge-provenance.md', folder / 'docs/merge-provenance.md')
    provenance = dict(product=PRODUCT, version=os.environ['PACKAGE_VERSION'], target=os.environ['NATIVE_TARGET'],
        head_sha=os.environ['GITHUB_SHA'], run_id=int(os.environ['GITHUB_RUN_ID']), run_attempt=int(os.environ['GITHUB_RUN_ATTEMPT']),
        native_host=os.environ['EXPECTED_HOST'], rustc=subprocess.check_output(['rustc', '-vV'], text=True),
        binary=identity(binary), lock_sha256=identity(ROOT / 'Cargo.lock')['sha256'],
        source_package=json.loads((delivery / 'package-verification.json').read_bytes()),
        native_dependencies=json.loads((delivery / 'native-dependencies.json').read_bytes()),
        accepted_delivery_base='f8c8f422ef8532701e3c70bbb019c72471b1a98d',
        historical_count_performance_run=40636288327, historical_ordinary_performance_run=41199697947,
        historical_ordinary_binary_sha256='926da17363eb37c469b075361c36930471c0a1425fc214ac881f3340cbbcac80',
        this_binary_remeasured=False, registry_publication=False, signed_release=False)
    save(folder / 'BUILD.json', provenance)
    for path in [folder, *folder.rglob('*')]:
        path.chmod(0o755 if path.is_dir() or path.name == PRODUCT else 0o644)
    archive = delivery / (stem + '.tar.gz')
    with tarfile.open(archive, 'w:gz') as output:
        output.add(folder, arcname=stem)
    source = Path(os.environ['SOURCE_PACKAGE'])
    (delivery / 'SHA256SUMS').write_text(''.join(identity(path)['sha256'] + '  ' + path.name + '\n' for path in (archive, source)))
    print(json.dumps(provenance, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('source', 'smoke', 'bundle'))
    parser.add_argument('--binary', type=Path)
    args = parser.parse_args()
    if args.action == 'source':
        source()
    elif args.action == 'smoke':
        smoke(args.binary)
    else:
        bundle(args.binary)


if __name__ == '__main__':
    main()
