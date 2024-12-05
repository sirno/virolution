//! # require_deferred_drop
//!
//! This crate provides a procedural macro to require deferred drop for a type that implements
//! `DeferredDrop`.

extern crate proc_macro;

use proc_macro::TokenStream;
use quote::quote;
use syn::{parse_macro_input, ItemFn};

#[proc_macro_attribute]
pub fn require_deferred_drop(_attr: TokenStream, item: TokenStream) -> TokenStream {
    // Parse the input tokens into a syntax tree
    let input_fn = parse_macro_input!(item as ItemFn);

    // Extract parts of the function
    let attrs = input_fn.attrs;
    let vis = input_fn.vis;
    let sig = input_fn.sig;
    let block = input_fn.block;

    // Generate the new function body
    let expanded = quote! {
        #(#attrs)*
        #vis #sig {
            self.require_deferred_drop();

            let __require_deferred_drop_result = (|| #block )();

            self.inquire_deferred_drop();

            __require_deferred_drop_result
        }
    };

    // Return the generated code as a TokenStream
    TokenStream::from(expanded)
}
