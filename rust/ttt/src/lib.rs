use std::{
    collections::HashMap,
    fmt::{
        Debug,
        Display,
        Formatter,
    },
    sync::mpsc::channel,
};

#[cfg(test)]
mod _test;

pub enum Error<'a> {
    Fail {
        test:   Option<Box<dyn Test + 'a>>,
        reason: String,
    },
}

impl<'a> Debug for Error<'a> {
    fn fmt(
        &self,
        f: &mut Formatter<'_>,
    ) -> std::fmt::Result {
        match self {
            Error::Fail {
                test,
                reason,
            } => {
                match test {
                    Some(t) => f.write_str(&format!("Some({})", t)),
                    None => f.write_str("None"),
                }?;
                write!(f, ": {reason}")
            }
        }
    }
}

pub trait Test: Display + Send {
    fn test<'t>(self: Box<Self>) -> Result<(), Error<'t>>;
}

pub struct Suite<'s> {
    name:  String,
    tests: HashMap<String, Box<dyn Test + 's>>,
}

impl<'s> Suite<'s> {
    pub fn new(name: String) -> Self {
        Self {
            name,
            tests: HashMap::new(),
        }
    }

    pub fn with_name(name: &str) -> Self {
        Self::new(name.to_string())
    }

    pub fn boxed(self) -> Box<Self> {
        Box::new(self)
    }

    fn insert(
        &mut self,
        test: Box<dyn Test + 's>,
    ) -> Option<Box<dyn Test + 's>> {
        self.tests.insert(test.to_string(), test)
    }

    pub fn append<T>(
        self,
        test: T,
    ) -> Self
    where
        T: Test + 's,
    {
        let mut suite = self;
        if suite.insert(Box::new(test)).is_some() {
            panic!("tests must have unique names")
        };
        suite
    }
}

impl<'s> Display for Suite<'s> {
    fn fmt(
        &self,
        f: &mut Formatter<'_>,
    ) -> std::fmt::Result {
        write!(f, "{}", self.name)
    }
}

impl<'s> Test for Suite<'s> {
    fn test<'t>(self: Box<Self>) -> Result<(), Error<'t>> {
        let tests = self.tests;
        let (tx, rx) = channel();
        rayon::scope(|s| {
            for (name, test) in tests {
                s.spawn(|_| {
                    tx.send((name, test.test())).expect("channel disconnected");
                });
            }
        });

        // Drop the last instance of tx end of the channel,
        // otherwise the iterator below will hang indefinitely.
        drop(tx);

        let suite_name = self.name;
        let mut failed = Suite::new(format!("failed from {suite_name}"));
        let mut all_good = Some(());
        for (name, res) in rx {
            match res {
                Ok(()) => eprintln!("{name} OK"),
                Err(Error::Fail {
                    reason,
                    test,
                }) => {
                    all_good = None;
                    eprintln!("{name} FAILED: {reason}");
                    if let Some(t) = test {
                        failed.tests.insert(t.to_string(), t);
                    }
                }
            }
        }
        all_good.ok_or(Error::Fail {
            test:   (!failed.tests.is_empty()).then_some(Box::new(failed)),
            reason: format!("failed in {suite_name}"),
        })
    }
}
