use std::{
    ffi::CString,
    path::Path,
};

mod ffi;

#[derive(Debug)]
enum Error {
    FileOpen,
}

struct DataHandle(ffi::dataid_t);

impl DataHandle {
    fn open(path: &Path) -> Result<Self, Error> {
        let path_str = path.to_str().ok_or(Error::FileOpen)?;
        let filename = CString::new(path_str).map_err(|_| Error::FileOpen)?;
        let file_id = unsafe { ffi::data_file_open(filename.as_ptr()) };
        if file_id == ffi::DATA_INVALID_OBJID {
            return Err(Error::FileOpen);
        }

        Ok(Self(file_id))
    }
}

impl Drop for DataHandle {
    fn drop(&mut self) {
        unsafe {
            ffi::data_file_close(self.0);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::path::{
        Path,
        PathBuf,
    };

    use super::*;

    #[test]
    fn data_file_open() {
        let data_dir = PathBuf::from("./dat");
        let filename = data_dir.join("./simul_H2_2.h5");
        let handle = DataHandle::open(&filename).unwrap();
    }
}
