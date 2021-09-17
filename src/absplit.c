/************************************************************************/
/**

   Program:    absplit
   \file       absplit.c
   
   \version    V0.1
   \date       17.09.21   
   \brief      Split an antibody PDB file into Fvs with antigens
   
   \copyright  (c) UCL / Prof. Andrew C. R. Martin 2021
   \author     Prof. Andrew C. R. Martin
   \par
               Institute of Structural & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   \par
               andrew@bioinf.org.uk
               andrew.martin@ucl.ac.uk
               
**************************************************************************

   This code is distributed under the conditions of the 
   GNU GENERAL PUBLIC LICENSE (GPL) Version 3

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

*************************************************************************/
/* Includes
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bioplib/SysDefs.h"

/************************************************************************/
/* Defines and macros
*/
#define MAXBUFF 240

/************************************************************************/
/* Globals
*/

/************************************************************************/
/* Prototypes
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile);


/************************************************************************/
int main(int argc, char **argv)
{
   char infile[MAXBUFF];
   
   if(ParseCmdLine(argc, argv, infile))
   {

   }
   
   
   return(0);
}


/************************************************************************/
/*>BOOL ParseCmdLine(int argc, char **argv, char *infile)
   ----------------------------------------------------------------------
*//**

   \param[in]      argc        Argument count
   \param[in]      **argv      Argument array
   \param[out]     *infile     Input filename (or blank string)
   \return                     Success

   Parse the command line

-  17.09.21 Original    By: ACRM
*/
BOOL ParseCmdLine(int argc, char **argv, char *infile)
{
   argc--;
   argv++;
   
   infile[0] = '\0';
   
   while(argc)
   {
      if(argv[0][0] == '-')
      {
         switch(argv[0][1])
         {
         case 'h':
            return(FALSE);
            break;
         default:
            return(FALSE);
            break;
         }
      }
      else
      {
         /* Check that there is only 1 argument left                    */
         if(argc != 1)
            return(FALSE);
         
         /* Copy the first to infile                                    */
         strcpy(infile, argv[0]);
         
         return(TRUE);
      }
      argc--;
      argv++;
   }
   
   return(TRUE);
}
