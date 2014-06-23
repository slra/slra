class Exception {
  static const size_t MSG_MAX = 200;

  char myMsg[MSG_MAX];
public:
  Exception( const char *msg, ... );
  const char *getMessage()  { return myMsg; }
};