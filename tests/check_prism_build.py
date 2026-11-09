"""Check PRISM profiling flags, schema selection, and request API isolation."""
import os
from pathlib import Path
import subprocess
import tempfile

root = Path(__file__).resolve().parents[1]
cc = os.environ.get('MPICC', 'mpicc')
cxx = os.environ.get('MPICXX', 'mpicxx')

with tempfile.TemporaryDirectory(prefix='qws-prism-build-') as directory:
  temp = Path(directory)
  harness = temp / 'schema.c'
  harness.write_text('#include <stdio.h>\nint rank=0;\n'
                     'FILE *para_outputfile=NULL;\n'
                     'void print_qws_timer_schema_(void);\n'
                     'int main(void) { print_qws_timer_schema_(); return 0; }\n')
  for name, flags, expected in (
      ('send', ['-DPROF_TARGET=TARGET_SEND'], ''),
      ('all', ['-DPROF_TARGET=TARGET_ALL'], ''),
      ('wait-only', ['-DPROF_TARGET=TARGET_OVERLAP_PROBE',
                     '-DOVERLAP_PROBE_FULL_TIMING', '-DOVERLAP_PROBE_WAIT_ONLY'],
       'probe_mode:wait_only'),
  ):
    executable = temp / name
    subprocess.run([cc, '-fopenmp', '-D_CHECK_TIMING', '-Werror=undef', *flags,
                    str(root/'timing.c'), str(harness), '-o', str(executable)],
                   check=True)
    output = subprocess.check_output([str(executable)], text=True)
    assert (expected in output if expected else output == ''), (name, output)
    if expected:
      assert 'hidden_recv_completion:not_observed' in output
    print('PASS schema', name, flush=True)

  common = [cxx, '-std=gnu++11', '-fopenmp', '-D_MPI_', '-D_CHECK_TIMING',
            '-D_NO_OMP_SINGLE', '-I'+str(root)]
  for name, flags in (
      ('send', ['-DPROF_TARGET=TARGET_SEND']),
      ('wait-only', ['-DPROF_TARGET=TARGET_OVERLAP_PROBE',
                     '-DOVERLAP_PROBE_FULL_TIMING', '-DOVERLAP_PROBE_WAIT_ONLY']),
  ):
    obj = temp / (name+'.o')
    subprocess.run(common+flags+['-c', str(root/'qws_xbound_mpi.cc'),
                                '-o', str(obj)], check=True)
    symbols = subprocess.check_output(['nm', '-u', str(obj)], text=True)
    for forbidden in ('MPI_Test', 'MPI_Cancel', 'MPI_Start', 'MPI_Wait',
                      'MPI_Send_init', 'MPI_Recv_init', 'MPI_Request_free'):
      assert not any(line.split()[-1] == forbidden for line in symbols.splitlines()), symbols
    assert 'OMPI_Prequest_free' in symbols
    assert 'OMPI_Wait' in symbols
    wanted = 'OMPI_Start_no_sync' if 'QWS_PRISM_START OMPI_Start_no_sync' in (root/'qws_prism.h').read_text() else 'OMPI_Start'
    assert any(line.split()[-1] == wanted for line in symbols.splitlines())
    print('PASS extension request isolation', name, wanted, flush=True)

  result = subprocess.run(common+['-DPROF_TARGET=TARGET_OVERLAP_PROBE',
                                  '-fsyntax-only', str(root/'qws_xbound_mpi.cc')],
                          capture_output=True, text=True)
  assert result.returncode != 0 and 'PRISM has no OMPI_Test' in result.stderr, result.stderr
  print('PASS unsupported polling probe rejected', flush=True)

  command = subprocess.check_output([
      'make', '-Bn', 'timing.o', 'fugaku_benchmark=', 'compiler=openmpi-gnu',
      'arch=skylake', 'mpi=1', 'omp=1', 'rdma=', 'powerapi=',
      'profiler=timing', 'target=send', 'CFLAGS=-O2 -fopenmp'], cwd=root, text=True)
  assert '-DPROF_TARGET=TARGET_SEND' in command and '-D_NO_OMP_SINGLE' in command, command
  print('PASS target and master-thread flags survive CFLAGS override', flush=True)
