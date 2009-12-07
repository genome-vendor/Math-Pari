#include "genpari.h"
#include <string.h>

#ifdef HAS_STRICMP
  #define strncasecmp strnicmp
#endif

#ifdef __cplusplus
extern "C" {
  /* int strncasecmp(const char *s1, const char *s2, int n); */
}
#endif

#ifdef READLINE2
  #if !defined(READLINE)
    #define READLINE
  #endif
#endif

#ifdef READLINE
  #ifdef __cplusplus 
    extern "C" {
  #endif 
  #include <readline/readline.h>
  #include <readline/history.h> 
  char **pari_completion (char *text, int start, int end);
  char *command_generator (char *text, int  state);
  void initialize_readline ();
  #ifdef __cplusplus 
  }
  #endif 
#endif

#ifdef _READLINE_H_

/*******************************************************************/
/*                                                                 */
/*                  Interface to Readline Completion               */
/*                                                                 */
/*******************************************************************/

/* Attempt to complete on the contents of TEXT.  START and END show the
   region of TEXT that contains the word to complete.  We can use the
   entire line in case we want to do some simple parsing.  Return the
   array of matches, or NULL if there aren't any. */

static entree* pe_compl;

#ifdef __cplusplus
extern "C" {
#endif

  typedef char** (*CMFUNC)(char*, char* (*)());

char **
pari_completion (char *text, int start, int end)
{
  char **matches;
  int first=0;
      /*  char *command_generator (); */

  matches = (char **)NULL;

  /* If the line does not begin in \, then it is a command
     to complete.  Otherwise it is the name of a file in the current
     directory. */

  while (rl_line_buffer[first] && isspace(rl_line_buffer[first])) first++;
  if (rl_line_buffer[first] != '\\')
  {
    matches=(char**)((CMFUNC)completion_matches)
      (text, (char* (*)())command_generator);
    
  /* Do some extra work to find out arguments if only one match found */
    if (matches && !matches[1] 
        && rl_line_buffer[end] != '(' && rl_line_buffer[first] != '?')
    {
      int n=0,len=strlen(pe_compl->name), addlen=0; 
      char *end=NULL, *beg;

      if (EpVALENCE(pe_compl)<100 && pe_compl>=fonctions
          && pe_compl<fonctions+NUMFUNC) n=pe_compl-fonctions;
      if(pe_compl->help && !strncmp(pe_compl->help, pe_compl->name, len)
         && pe_compl->help[len] == '(' 
	 && (end = strchr(pe_compl->help, ')')))
      {
	beg = pe_compl->help + len;
	addlen = end - beg + 1;
      }
      else if (EpVALENCE(pe_compl) != 200) {beg = "()"; addlen = 2;}
      if (addlen)
      {
        matches[0] = (char*)realloc(matches[0], 
				    strlen(matches[0]) + 1 + addlen);
        strncat(matches[0], beg, addlen);
      }
    }
  }
  return (matches);
}
#ifdef __cplusplus
}
#endif


/* Generator function for command completion.  STATE lets us know whether
   to start from scratch; without any state (i.e. STATE == 0), then we
   start at the top of the list. */
char *
command_generator (char *text, int  state)
{
  static int n, len, junk;
  static entree* ep;
  char *name;

  /* If this is a new word to complete, initialize now.  This includes
     saving the length of TEXT for efficiency, and initializing the index
     variable to 0. Since file completion and symbol completion use
     different word boundaries, put a new boundary into junk. */
  if (!state)
  {
    n = 0; ep=hashtable[n];
    len = strlen (text);
    junk=len-1;
    while (junk >= 0 && isalnum(text[junk])) junk--;  
    junk++;
  }

  /* Return the next name which partially matches from the command list. */
  while (!(ep && !strncasecmp(ep->name,text+junk,len-junk)) && (ep || ++n<TBLSZ))
  {
    ep? (ep=ep->next): (ep=hashtable[n]);
  }
  if (ep)
  {
    name=strncpy((char*)malloc(strlen(ep->name)+1+junk),text,junk);
    strcpy(name+junk,ep->name); pe_compl=ep;
    ep=ep->next;
    return(name);
  }

  /* If no names matched, then return NULL. */
  return ((char *)NULL);
}

#ifdef READLINE2

int
rl_short_help(int count, int key)
{
    int off=rl_point, n, off1, found=0, point=off, len;
    entree *ep;
    char *u, *v;

    while (off && isalnum(rl_line_buffer[off-1])) {
	off--;
    }
    for (off1=off, n=0; isalnum(rl_line_buffer[off1]); off1++)
	n = n << 1 ^ rl_line_buffer[off1];
    len = off1 - off;
    if (n < 0) n = -n; n %= TBLSZ;
    for(ep = hashtable[n]; ep; ep = ep->next)
	{
	    if (len == strlen(ep->name) 
		&& !strncasecmp(ep->name, rl_line_buffer + off, len)) {
		found=1;
		break;
	    }
	}
    if (!found || !ep->help) { 
	ding(); 
/*      rl_message("Unknown name `%*s', press a key! ",
        off1-off,rl_line_buffer[off]); 
        rl_get_char();
        rl_clear_message(); */
	return 0;
    }
    rl_point=0;
    ((void(*)(const char*))rl_message)("");
    if (count>0) {
	fprintf(rl_outstream, "%s.\n", ep->help);
    } else {
	char command[80];
	sprintf(command, "gphelp \"%.70s\"", ep->name);
	system(command) && ding();
    }
    /*fflush(rl_outstream);*/
    rl_point=point;
    rl_on_new_line();
    rl_clear_message(); 
    return 0;
}

int rl_long_help(int count, int key)
{
  return rl_short_help(-1,key);
}

#endif /* defined(READLINE2) */

/* Tell the GNU Readline library how to complete.  We want to try to complete
   on command names if the first char in the line is not \\, or on filenames
   if it is. */

void
initialize_readline ()
{
      /*  char **pari_completion (); */

      /* Allow conditional parsing of the ~/.inputrc file. */
  rl_readline_name = "Pari-GP";

      /* Tell the completer that we want a crack first. */
#ifdef READLINE2
  rl_attempted_completion_function = (CPPFunction *)pari_completion;
  ((void(*)(const char*,Function*,int))rl_add_defun)
    ("short-help",(Function*)rl_short_help,-1);
  ((void(*)(const char*,Function*,int))rl_add_defun)
    ("long-help",(Function*)rl_long_help,-1);
  ((void(*)(int,char*,char*,Keymap))rl_generic_bind)
    (ISFUNC,"h",(char*)rl_short_help,emacs_meta_keymap);
  ((void(*)(int,char*,char*,Keymap))rl_generic_bind)
    (ISFUNC,"H",(char*)rl_long_help,emacs_meta_keymap);
  ((void(*)(int,char*,char*,Keymap))rl_generic_bind)
    (ISFUNC,"h",(char*)rl_short_help,vi_movement_keymap);
  ((void(*)(int,char*,char*,Keymap))rl_generic_bind)
    (ISFUNC,"H",(char*)rl_long_help,vi_movement_keymap);
#ifdef EMACS_DOS_KEYMAP
  /* F1 */
  ((void(*)(int,char*,char*,Keymap))rl_generic_bind)
    (ISFUNC,";",(char*)rl_short_help,emacs_dos_keymap);
  /* Shift-F1 */
  ((void(*)(int,char*,char*,Keymap))rl_generic_bind)
    (ISFUNC,"T",(char*)rl_long_help,emacs_dos_keymap);
#endif
#else
  rl_attempted_completion_function = (Function *)pari_completion;
#endif
}

#else
  int aBogusVariable=1;
#endif /* _READLINE_H_ */
