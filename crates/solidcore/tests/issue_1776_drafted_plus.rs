//! tangency-labs/solidcore#1776 (WR-587): `SmCurve::ComputeCrvsNormal` sums the cross
//! products of consecutive chords, so each turn counts by its chords' lengths. A plus sign
//! of arms 3 long and 1 wide, drawn counter-clockwise, has 8 convex turns between chords
//! of 1 and 3 (+24) and 4 reflex turns between chords of 3 and 3 (-36): its normal points
//! down, and `CreateTaperExtrude`, offsetting about that normal, drafts the plus against
//! the requested direction.
//!
//! With `Kernel::fixes.crvs_normal_area_vector` the normal is the area vector of the
//! polygon through the samples (Newell's method), which follows the loop's orientation:
//! up for the counter-clockwise plus. Without it, the pinned behaviour stays: outward and
//! inward drafts swap their volumes.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing
)] // test helpers
use solidcore::curve::curve::p3::compute_crvs_normal;
use solidcore::curve::curve::Curve;
use solidcore::curve::line::Line;
use solidcore::geom::attribute::{AttributeStore, ObjectKey};
use solidcore::kernel::{Kernel, KernelOptions};
use solidcore::mass::integrator::properties::{brep_volume, PropertyAttributes};
use solidcore::sweep::api_brep::{
    create_taper_extrude_with_corner, OffsetCornerType, TaperExtrude,
};
use solidcore::topology::brep::p2::is_manifold_solid;
use solidcore::Vec3;

/// The height of the drafts: an inward offset of 2 tan 10 = 0.35 keeps the arms (1 wide)
/// open to the far end.
const HEIGHT: f64 = 2.0;

/// The plus sign, counter-clockwise (area 13, perimeter 28).
fn plus() -> Vec<Curve> {
    let pts = [
        (0.5, -3.5),
        (0.5, -0.5),
        (3.5, -0.5),
        (3.5, 0.5),
        (0.5, 0.5),
        (0.5, 3.5),
        (-0.5, 3.5),
        (-0.5, 0.5),
        (-3.5, 0.5),
        (-3.5, -0.5),
        (-0.5, -0.5),
        (-0.5, -3.5),
    ]
    .map(|(x, y)| Vec3::new(x, y, 0.0));
    (0..pts.len())
        .map(|i| {
            Curve::from(
                Line::create_line_segment(3, pts[i], pts[(i + 1) % pts.len()], None).unwrap(),
            )
        })
        .collect()
}

fn kernel(fix: bool) -> Kernel {
    let mut k = Kernel::new(&KernelOptions::default());
    k.fixes.crvs_normal_area_vector = fix;
    k
}

/// The plus drafted 10 degrees (positive inward, negative outward), mitred corners: a
/// closed manifold solid's volume.
fn drafted(k: &mut Kernel, draft_deg: f64) -> f64 {
    let curves = plus();
    let refs: Vec<Option<&Curve>> = curves.iter().map(Some).collect();
    let result = create_taper_extrude_with_corner(
        k,
        &refs,
        HEIGHT,
        draft_deg,
        3,
        OffsetCornerType::LinearExtension,
    )
    .unwrap();
    let TaperExtrude::Solid(mut b) = result else {
        panic!("refused: {result:?}");
    };
    assert!(is_manifold_solid(&b).unwrap());
    let store = AttributeStore::default();
    let key_of = |_| ObjectKey(0);
    let attributes = PropertyAttributes {
        store: &store,
        key_of: &key_of,
    };
    brep_volume(k, &mut b, 1e-6, &attributes).unwrap()
}

/// A rectilinear profile of area a and perimeter p, with c convex and r reflex corners,
/// drafted by offset o over the height h with mitred walls:
/// h (a + s p o / 2 + (c - r) o^2 / 3), s = +1 outward, -1 inward.
fn plus_volume(outward: bool) -> f64 {
    let o = HEIGHT * 10.0_f64.to_radians().tan();
    let s = if outward { 1.0 } else { -1.0 };
    HEIGHT * (13.0 + s * 28.0 * o / 2.0 + 4.0 * o * o / 3.0)
}

fn assert_rel(a: f64, b: f64) {
    assert!((a - b).abs() <= 1e-9 * a.abs().max(b.abs()), "{a} vs {b}");
}

#[test]
fn the_area_vector_fix_is_off_by_default() {
    assert!(!Kernel::default().fixes.crvs_normal_area_vector);
}

#[test]
fn the_plus_normal_points_down_by_default_and_up_with_the_fix() {
    let curves = plus();
    assert!(compute_crvs_normal(&kernel(false), &curves).unwrap().z < -0.99);
    assert!(compute_crvs_normal(&kernel(true), &curves).unwrap().z > 0.99);
}

#[test]
fn a_convex_loop_s_normal_is_unchanged_by_the_fix() {
    let pts = [(0.0, 0.0), (4.0, 0.0), (4.0, 1.0), (0.0, 1.0)].map(|(x, y)| Vec3::new(x, y, 0.0));
    let curves: Vec<Curve> = (0..4)
        .map(|i| Curve::from(Line::create_line_segment(3, pts[i], pts[(i + 1) % 4], None).unwrap()))
        .collect();
    let before = compute_crvs_normal(&kernel(false), &curves).unwrap();
    let after = compute_crvs_normal(&kernel(true), &curves).unwrap();
    assert!(before.z > 0.99 && after.z > 0.99);
}

#[test]
fn with_the_fix_the_drafts_lean_as_asked() {
    let mut k = kernel(true);
    assert_rel(drafted(&mut k, -10.0), plus_volume(true));
    assert_rel(drafted(&mut k, 10.0), plus_volume(false));
}

#[test]
fn without_the_fix_outward_and_inward_swap() {
    let mut k = kernel(false);
    assert_rel(drafted(&mut k, -10.0), plus_volume(false));
    assert_rel(drafted(&mut k, 10.0), plus_volume(true));
}
