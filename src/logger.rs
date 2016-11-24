extern crate log;
extern crate colored;

use log::{LogRecord, LogLevel, LogMetadata, SetLoggerError, LogLevelFilter};
use colored::*;

pub struct Logger {
    level: LogLevelFilter,
}

impl log::Log for Logger {
    fn enabled(&self, metadata: &LogMetadata) -> bool {
        metadata.level() <= LogLevel::Trace
    }

    fn log(&self, record: &LogRecord) {
        if record.level() <= self.level {
            println!("{}",
                     match record.level() {
                         LogLevel::Error => { format!("{} {}", "✖".red(),    record.args().to_string().red().bold()) }
                         LogLevel::Warn  => { format!("{} {}", "⚠".yellow(), record.args().to_string().yellow()) }
                         LogLevel::Info  => { format!("{} {}", "⬨".cyan(),   record.args()) }
                         LogLevel::Trace => { format!("{} {}", "✓".green(),  record.args().to_string().green()) }
                         LogLevel::Debug => { format!("{} {}", "❖".blue(),   record.args().to_string().blue()) }
                     });
        }
    }
}

impl Logger {
    pub fn init(level: LogLevelFilter) -> Result<(), SetLoggerError> {
        log::set_logger(|max_log_level| {
            max_log_level.set(level);
            Box::new(Logger {level: level})
        })
    }
}
