extern crate colored;
extern crate log;

use self::log::{Record, Level, Metadata, SetLoggerError, LevelFilter};
use self::colored::*;

pub struct Logger {
}

impl log::Log for Logger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() >= Level::Trace
    }

    fn log(&self, record: &Record) {
        println!("{}",
                 match record.level() {
                     Level::Error => {
                         format!("{} {}", "✖".red(), record.args().to_string().red().bold())
                     }
                     Level::Warn => {
                         format!("{} {}", "⚠".yellow(), record.args().to_string().yellow())
                     }
                     Level::Info => {
                         format!("{}", record.args())
                     }
                     Level::Trace => {
                         format!("{} {}", "▷".cyan(), record.args())
                     }
                     Level::Debug => {
                         format!("{} {}", "❖".blue(), record.args().to_string().blue())
                     }
                 });
    }

    fn flush(&self) {}
}

impl Logger {
    pub fn init(level: LevelFilter) -> Result<(), SetLoggerError> {
        log::set_boxed_logger(Box::new(Logger{})).map(|()| log::set_max_level(level))
    }
}
