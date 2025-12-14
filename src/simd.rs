use std::sync::LazyLock;

use pulp::Arch;

pub(crate) static ARCH: LazyLock<Arch> = LazyLock::new(Arch::new);

pub(crate) fn simd_sum(values: &[f32]) -> f32 {
    let mut total = 0f32;
    ARCH.dispatch(|| {
        for x in values {
            total += x;
        }
    });
    total
}
