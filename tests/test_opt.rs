// [[file:../optim.note::f15831bf][f15831bf]]
use gosh_core::*;
use gut::prelude::*;

#[test]
fn test_opt() -> Result<()> {
    use gchemol::prelude::*;
    use gchemol::Molecule;
    use gosh_model::LennardJones;
    use gosh_optim::{optimize_geometry_iter, Optimizer};

    let filename = "tests/files/LennardJones/LJ38r.xyz";
    let mut mol = Molecule::from_file(filename)?;
    let mut lj = LennardJones::default();
    lj.derivative_order = 1;

    let _ = Optimizer::default().optimize_geometry(&mut mol, &mut lj)?;

    // iterator interface
    let steps = optimize_geometry_iter(&mut mol, &mut lj);
    for p in steps.take(10) {
        dbg!(p.fmax, p.ncalls);
    }

    Ok(())
}
// f15831bf ends here
