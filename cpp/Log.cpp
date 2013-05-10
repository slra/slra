#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <cstdarg>

#include "slra.h"

void Log::lprintf( Level level, char *format, ... ) {
  va_list vl;
  va_start(vl, format);  
  Log * log = getLog();

  if (level  <= log->myMaxLevel) { 
    log->myMsg[MSG_MAX-1] = 0;
    vsnprintf(log->myMsg, MSG_MAX-1, format, vl); 
  
    PRINTF(log->myMsg);
  }
}

void Log::lprintf( char *format, ... ) {
  va_list vl;
  va_start(vl, format);  
  Log * log = getLog();

  log->myMsg[MSG_MAX-1] = 0;
  vsnprintf(log->myMsg, MSG_MAX-1, format, vl); 
  PRINTF(log->myMsg);
  FLUSH();
}

  
void Log::setMaxLevel( Level maxLevel ) {
  getLog()->myMaxLevel = maxLevel;
}
  
void Log::str2DispLevel( const char *str )  {
  char *str_disp[] = { "off", "final", "notify", "iter" };
  size_t i;
  
  
  for (i = 0; i < sizeof(str_disp) / sizeof(str_disp[0]); i++) {
    if (strcmp(str_disp[i], str) == 0) {
      getLog()->myMaxLevel = (Level)i;
    }
  }
    
  if (i < 0) {
    MYWARNING("Ignoring optimization option 'disp'. Unrecognized value.\n");
  }
}

Log *Log::myLogInstance = NULL;
  
Log::Level Log::getMaxLevel() {
  return getLog()->myMaxLevel;
}

Log *Log::getLog() {
  if (myLogInstance == NULL) {
    myLogInstance = new Log;
  }
  
  return myLogInstance;
}

void Log::deleteLog() {
  if (myLogInstance != NULL) {
    delete myLogInstance;
    myLogInstance = NULL;
  }
}
