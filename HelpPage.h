#ifndef HELPPAGE_H
#define HELPPAGE_H

// -h or -v
int HelpPage(FILE *fp);

void ErrMsg(double messageno);

// general warnings
void WarMsg(double messageno);

float myatof(char *str);

void Sizebreak(char *str);

#endif