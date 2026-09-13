"""Graph-based node-shape classification for topology compatibility."""
from __future__ import annotations
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from .model import MonomerSpec
SHAPE_SQUARE='square'; SHAPE_RECTANGULAR='rectangular'; SHAPE_TETRAHEDRAL='tetrahedral'; SHAPE_UNKNOWN='unknown'; SHAPE_C4_TD='c4_td'; SHAPE_D2H='d2h'
KNOWN_SHAPE_LABELS=frozenset({SHAPE_SQUARE,SHAPE_RECTANGULAR,SHAPE_TETRAHEDRAL})
@dataclass(frozen=True)
class NodeShapeSignature: label: str
def unknown_signature(): return NodeShapeSignature(SHAPE_UNKNOWN)
def classify_monomer_graph_shape(spec: MonomerSpec) -> str:
    if len(spec.motifs)!=4 or not spec.bonds: return SHAPE_UNKNOWN
    adj={}
    for a,b,_ in spec.bonds: adj.setdefault(a,set()).add(b); adj.setdefault(b,set()).add(a)
    sites=[]
    for motif in spec.motifs:
        members=set(motif.atom_ids); candidates=[a for a in members if a<len(spec.atom_symbols) and spec.atom_symbols[a]!='H']
        if not candidates:return SHAPE_UNKNOWN
        sites.append(max(candidates,key=lambda a:sum(n not in members for n in adj.get(a,()))))
    if len(set(sites))!=4:return SHAPE_UNKNOWN
    distances=[]
    for i,source in enumerate(sites):
        seen={source:0}; queue=[source]
        for node in queue:
            for n in adj.get(node,()):
                if n not in seen: seen[n]=seen[node]+1; queue.append(n)
        for target in sites[i+1:]:
            if target not in seen:return SHAPE_UNKNOWN
            distances.append(seen[target])
    colors={n:(spec.atom_symbols[n],len(adj.get(n,()))) for n in adj}
    for _ in range(max(1,len(adj))):
        keys={n:(colors[n],tuple(sorted(colors.get(x,()) for x in adj.get(n,())))) for n in adj}
        palette={k:i for i,k in enumerate(sorted(set(keys.values()),key=repr))}; colors={n:palette[k] for n,k in keys.items()}
    orbit_counts=sorted(Counter(colors.get(s,-1) for s in sites).values()); counts=sorted(Counter(distances).values())
    if counts==[6] or (orbit_counts==[4] and len(set(distances))<=2): return SHAPE_C4_TD
    if counts in ([2,4],[2,2,2]): return SHAPE_D2H
    return SHAPE_UNKNOWN
def classify_monomer_node_shape(spec):
    family=classify_monomer_graph_shape(spec)
    return NodeShapeSignature(SHAPE_TETRAHEDRAL if family==SHAPE_C4_TD else SHAPE_RECTANGULAR if family==SHAPE_D2H else SHAPE_UNKNOWN)
@lru_cache(maxsize=None)
def classify_topology_node_shape(topology_id): return NodeShapeSignature({'sql':SHAPE_SQUARE,'kgm':SHAPE_RECTANGULAR,'dia':SHAPE_TETRAHEDRAL}.get(topology_id,SHAPE_UNKNOWN))
def shapes_compatible(a,b):
    if a.label not in KNOWN_SHAPE_LABELS or b.label not in KNOWN_SHAPE_LABELS:return None
    return a.label==b.label
