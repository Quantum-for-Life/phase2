mod circ;

mod quest_sys;

mod linen;
mod rayon;

/// The impl of `SendPtr` was taken from [`rayon`] crate.
/// Rayon is distributed under MIT License.
///
/// We need to transmit raw pointers across threads. It is possible to do this
/// without any unsafe code by converting pointers to usize or to `AtomicPtr`<T>
/// then back to a raw pointer for use. We prefer this approach because code
/// that uses this type is more explicit.
///
/// Unsafe code is still required to dereference the pointer, so this type is
/// not unsound on its own, although it does partly lift the unconditional
/// !Send and !Sync on raw pointers. As always, dereference with care.
///
/// [`rayon`]: https://crates.io/crates/rayon
#[derive(Debug)]
#[repr(transparent)]
pub(crate) struct SendPtr<T>(*mut T);

// SAFETY: !Send for raw pointers is not for safety, just as a lint
unsafe impl<T: Send> Send for SendPtr<T> {}

// SAFETY: !Sync for raw pointers is not for safety, just as a lint
unsafe impl<T: Send> Sync for SendPtr<T> {}

impl<T> SendPtr<T> {
    // Helper to avoid disjoint captures of `send_ptr.0`
    #[allow(dead_code)]
    #[must_use]
    pub fn get(self) -> *mut T {
        self.0
    }
}

// Implement Copy without the T: Copy bound from the derive
impl<T> Copy for SendPtr<T> {}

// Implement Clone without the T: Clone bound from the derive
impl<T> Clone for SendPtr<T> {
    fn clone(&self) -> Self {
        *self
    }
}
