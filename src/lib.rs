// [[file:../optim.note::a74a585a][a74a585a]]
use gosh_core::*;
use vecfx::*;

use gut::prelude::*;
// a74a585a ends here

// [[file:../optim.note::2e984082][2e984082]]
mod opt;
mod optimization;
mod potential;
mod vars;
// 2e984082 ends here

// [[file:../optim.note::33bebce4][33bebce4]]
pub use opt::*;
pub use potential::{Dynamics, EvaluatePotential, PotentialOutput};
// 33bebce4 ends here

// [[file:../optim.note::242ad86a][242ad86a]]
#[cfg(feature = "adhoc")]
/// Docs for local mods
pub mod docs {
    macro_rules! export_doc {
        ($l:ident) => {
            pub mod $l {
                pub use crate::$l::*;
            }
        };
    }

    export_doc!(potential);
    export_doc!(opt);
    export_doc!(vars);
}
// 242ad86a ends here
