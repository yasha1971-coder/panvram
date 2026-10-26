import sys, re, torch
torch.cuda.current_device = lambda: "cpu"   # test only: the model hardcodes the current CUDA device
sys.path.insert(0,'.')
EXTRA = sys.argv[1:]
src=open('/tmp/claude-1002/m3b/m3b.py').read()
common=re.search(r'COMMON = f"""(.*?)"""', src, re.S).group(1)
common=common.replace("{SEQ}","256").replace("{MBS}","2").replace("{GBS}","2").replace("{ITERS}","200")
sys.argv=['pretrain_gpt.py']+[x for x in common.split() if x!='--bf16']+['--data-path','panvram:/content/cohort','--use-cpu-initialization']+EXTRA
from megatron.training.arguments import parse_args, validate_args
from megatron.training import global_vars
a=parse_args(); validate_args(a,{}); global_vars.set_args(a)
torch.distributed.init_process_group("gloo", world_size=1, rank=0)
from megatron.core import parallel_state
parallel_state.initialize_model_parallel(1,1)
import megatron.core.transformer.utils as tu
def _mask(sq, sk=None):
    sk = sq if sk is None else sk
    return torch.triu(torch.ones(sq, sk, dtype=torch.bool), diagonal=1)
tu.get_default_causal_mask = _mask
import megatron.core.fusions.fused_softmax as fs
fs.get_default_causal_mask = _mask
import contextlib, megatron.core.tensor_parallel.random as tr
class _T:
    def fork(self, *a, **k): return contextlib.nullcontext()
tr.get_cuda_rng_tracker = lambda *a, **k: _T()
import megatron.core.tensor_parallel as tp_; tp_.get_cuda_rng_tracker = tr.get_cuda_rng_tracker
import gpt_builders, model_provider
from megatron.training.tokenizer import build_tokenizer
global_vars._GLOBAL_TOKENIZER = build_tokenizer(a)
m = model_provider.model_provider(gpt_builders.gpt_builder, pre_process=True, post_process=True)
print('MODEL BUILT', sum(p.numel() for p in m.parameters()))
tok = torch.randint(5,10,(2,256)); pos = torch.arange(256).expand(2,256)
out = m(tok, pos, None, labels=tok)
out.float().mean().backward()
print('FORWARD+BACKWARD OK', tuple(out.shape), float(out.float().mean()))
