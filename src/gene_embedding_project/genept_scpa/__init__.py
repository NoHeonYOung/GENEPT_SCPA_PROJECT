"""GenePT × SCPA experiment package.

Only protocol/configuration utilities are active in Phase 0. Analysis modules are
reserved so later phases do not mix with unrelated pathway-recovery code.
"""

from .config import ProtocolConfig, load_config

__all__ = ["ProtocolConfig", "load_config"]

