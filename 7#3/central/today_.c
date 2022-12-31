/*
*   C SUBROUTINE: today_
*
*   PURPOSE: Used by fortran xxtoday subroutine.
*
*   NOTES:
*
*   AUTHOR:  Martin O'Mullane
*   DATE:    19-04-2005
*
*   UPDATE:
*
*/

#include <time.h>

long today_(void)
{
  int year, month, day;
  long value;

  time_t tnow;
  struct tm *tmnow;

  time(&tnow);
  tmnow = localtime(&tnow);

  year  = tmnow->tm_year + 1900;
  month = tmnow->tm_mon + 1;
  day   = tmnow->tm_mday;

  value = day + month*100 + year*10000;

  return value;

}
