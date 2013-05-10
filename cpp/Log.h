/** Class for logging */
class Log {
public:  
  enum Level {
    LOG_LEVEL_OFF = 0,
    LOG_LEVEL_FINAL,
    LOG_LEVEL_NOTIFY,
    LOG_LEVEL_ITER
  };

  static void lprintf( char *format, ... );  
  static void lprintf( Level level, char *format, ... );  

  /** Initialize disp field from string */
  static void str2DispLevel( const char *str );
  static void setMaxLevel( Level maxLevel );
  static Level getMaxLevel();

  static void deleteLog();
  
private:
  static const size_t MSG_MAX = 200;
  char myMsg[MSG_MAX];
  Level myMaxLevel;

  static Log *getLog();


  Log() : myMaxLevel(LOG_LEVEL_NOTIFY) {};
  ~Log() {};
  static Log *myLogInstance;
};
