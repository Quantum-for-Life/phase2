use std::fmt::{
    Display,
    Formatter,
};

use crate::{
    Error,
    Suite,
    Test,
};

fn add(
    x: usize,
    y: usize,
) -> usize {
    x + y
}

#[derive(Debug)]
struct TestAdd(usize, usize);

impl Display for TestAdd {
    fn fmt(
        &self,
        f: &mut Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Test for TestAdd {
    fn test<'t>(self: Box<Self>) -> Result<(), Error<'t>> {
        if add(self.0, self.1) == self.0 + self.1 {
            Ok(())
        } else {
            Err(Error::Fail {
                test:   Some(self),
                reason: "adding error".to_string(),
            })
        }
    }
}

#[derive(Debug)]
struct TestFail(u8);

impl Display for TestFail {
    fn fmt(
        &self,
        f: &mut Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Test for TestFail {
    fn test<'t>(mut self: Box<Self>) -> Result<(), Error<'t>> {
        Err(Error::Fail {
            test:   if self.0 > 0 {
                self.0 -= 1;
                Some(self)
            } else {
                None
            },
            reason: "always fails".to_string(),
        })
    }
}

#[test]
fn test_suite() {
    let suite_1 = Suite::with_name("1")
        .append(TestAdd(0, 0))
        .append(TestAdd(1, 1))
        .append(TestFail(1))
        .append(TestAdd(2, 2));
    let suite_2 = Suite::with_name("2")
        .append(TestAdd(0, 0))
        .append(TestAdd(1, 1))
        .append(TestFail(1))
        .append(TestAdd(2, 2));
    let suite = Suite::with_name("3")
        .append(suite_1)
        .append(suite_2)
        .boxed();

    let failed = suite.test().unwrap_err();

    println!("{failed:?}")
}
