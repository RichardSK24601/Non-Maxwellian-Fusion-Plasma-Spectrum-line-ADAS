/*
* 
*   C SUBROUTINE: pers_f
* 
*   PURPOSE: Given a username returns the users real name from his
*            password entry.       
*
*   NOTES:  This routine is called by the fortran routine xxname
*
*   AUTHOR:  Martin O'Mullane (probably)
*
*   DATE:    11-08-98 (maybe)
*
*   VERSION :  1.1
*   DATE    :  11-08-98 (maybe)
*   UPDATE  :  Martin O'Mullane (probably)
*               - First release 
*
*   VERSION :  1.2
*   DATE    :  19-06-2006
*   UPDATE  :  Allan Whiteford
*               - Added missing documentation (guess at version history)
*               - Blank string before call to getpwnam
*               - Return early if getpwname fails
*
*/

#include <pwd.h>
void pers_f_ (char *ret, int retl, char *str, int slen)
{
  int i;
  int len;
  char *tmp = 0;
  struct passwd *ent;
  
  for (i=0; i<retl; i++) 
      *(ret+i) = ' '; 

  if (ent = getpwnam(str)) {
    tmp = ent->pw_gecos;
  }  else return;
  
  while (*ret++ = *tmp++) ;
  
  return;
}
