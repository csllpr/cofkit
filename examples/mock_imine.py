from cofkit import COFEngine, COFProject, Frame, MonomerSpec, ReactiveMotif


def build_project() -> COFProject:
    tri_amine = MonomerSpec(
        id="tapb",
        name="TAPB-like triamine",
        motifs=(
            ReactiveMotif(id="n1", kind="amine", atom_ids=(1,), frame=Frame.xy()),
            ReactiveMotif(id="n2", kind="amine", atom_ids=(2,), frame=Frame.yz()),
            ReactiveMotif(id="n3", kind="amine", atom_ids=(3,), frame=Frame.zx()),
        ),
    )
    tri_aldehyde = MonomerSpec(
        id="tfp",
        name="TFP-like trialdehyde",
        motifs=(
            ReactiveMotif(id="c1", kind="aldehyde", atom_ids=(4,), frame=Frame.xy()),
            ReactiveMotif(id="c2", kind="aldehyde", atom_ids=(5,), frame=Frame.yz()),
            ReactiveMotif(id="c3", kind="aldehyde", atom_ids=(6,), frame=Frame.zx()),
        ),
    )
    return COFProject(
        monomers=(tri_amine, tri_aldehyde),
        allowed_reactions=("imine_bridge",),
        target_dimensionality="2D",
        target_topologies=("hcb",),
    )


if __name__ == "__main__":
    engine = COFEngine()
    results = engine.run(build_project())
    best = results.top(1)[0]
    print("candidate:", best.id)
    print("score:", best.score)
    print("flags:", ", ".join(best.flags) or "<none>")
    print("events:", len(best.events))
    print("summary:", best.metadata["graph_summary"])
