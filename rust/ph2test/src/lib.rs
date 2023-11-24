use phase2::circ::CircEnv;

pub fn hello() {
    println!("Hello from ph2test!");

    let env = CircEnv::try_new().unwrap();
    env.report();
}
