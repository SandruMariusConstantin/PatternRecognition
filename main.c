#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct{
    unsigned char R;
    unsigned char G;
    unsigned char B;
}pixel;

typedef struct{
    unsigned int i;
    unsigned int j;
}fereastra;

typedef struct{
    double corelatiee;
    fereastra fe;
    unsigned int exista;
    pixel culoare;
}detectie;

unsigned int * XORSHIFT32(unsigned int seed, unsigned int n)
{

    unsigned int k, r, *R;
    R = (unsigned int) calloc (2*n, sizeof(unsigned int));
    r = R[0] = seed;

    for(k = 1; k < 2*n; k++)
    {
        r ^= r << 13;
        r ^= r >> 17;
        r ^= r << 5;
        R[k] = r;
    }

    return R;
}

pixel* incarca_imagine (char *nume_fisier_sursa)
{
   FILE *fin;
   unsigned int latime_img, inaltime_img;
   pixel *L;

   fin = fopen(nume_fisier_sursa, "rb");
   if(fin == NULL)
   {
       printf("Nu am putut gasi imaginea\n");
       return;
   }

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img,   sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);


   unsigned int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;



    fseek(fin, 54, SEEK_SET);
    unsigned int i,j;

    L = (pixel *) calloc (latime_img*inaltime_img, 3);
    if(L == NULL)
    {
        printf("Eroare la alocare");
        return;
    }

	for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++){
            fread(&(L[(inaltime_img - i - 1) * latime_img + j].B), 1, 1, fin);
            fread(&(L[(inaltime_img - i - 1) * latime_img + j].G), 1, 1, fin);
            fread(&(L[(inaltime_img - i - 1) * latime_img + j].R), 1, 1, fin);} ///matricea e un vector unde: index_vector=index_linie_matrice*nr_coloane+index_coloane_matrice
        fseek(fin, padding, SEEK_CUR);
    }

    fclose(fin);
    return L;
}

void salvare_imagine(char *nume_fisier_sursa, char *nume_fisier_destinatie, pixel *C)
{
    FILE *fin, *fout;
    unsigned int latime_img, inaltime_img, i, j;

    fin = fopen(nume_fisier_sursa, "rb");
    fout = fopen(nume_fisier_destinatie, "wb");

    if(fin == NULL)
    {
        printf("Nu am putut gasi imaginea\n");
        return 0;
    }

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);

   unsigned int padding;
   unsigned char padd = 0x00;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    fseek(fin,0,SEEK_SET);

	unsigned char c;
    for(i = 0; i < 54; i++)
    {
        fread(&c, 1, 1, fin);
		fwrite(&c, 1, 1, fout);
		fflush(fout);
	}
	fclose(fin);

	for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            fwrite(&(C[(inaltime_img - i - 1) * latime_img + j].B), 1, 1, fout); fflush(fout);
            fwrite(&(C[(inaltime_img - i - 1) * latime_img + j].G), 1, 1, fout); fflush(fout);
            fwrite(&(C[(inaltime_img - i - 1) * latime_img + j].R), 1, 1, fout); fflush(fout);
        }
        for(int pd = 0; pd < padding; pd++){
        fwrite(&padd, 1, 1, fout); fflush(fout);}
    }

    fclose(fout);
    printf("Imaginea a fost salvata\n\n");
}

void criptare(char *nume_fisier_sursa, char *nume_fisier_destinatie, char *cheie_secreta)
{
   ///Luam cheia
   unsigned int R0, SV;
   FILE *key;
   key = fopen(cheie_secreta, "r");
   if(key == NULL)
   {
       printf("Fisierul unde se afla cheia nu a putut fi deschis/n");
       return;
   }
   fscanf(key, "%u", &R0); fscanf(key, "%u", &SV);
   fclose(key);


   FILE *fin;
   fin = fopen(nume_fisier_sursa, "rb");
   unsigned int latime_img, inaltime_img, n;

   if(fin == NULL)
   {
       printf("Nu am putut gasi imaginea");
       return;
   }

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img,   sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);
   fclose(fin);

   n = latime_img * inaltime_img;

   pixel *L, *Lprim;
   unsigned int *R, i;

   R = XORSHIFT32 (R0, n);
   L = incarca_imagine (nume_fisier_sursa);

   Lprim = (pixel *) calloc (n, 3);
   unsigned int *perm = calloc (n, sizeof(unsigned int));

   ///PERMUTARE
   for(i = 0; i < n; i++) perm[i] = i;
   unsigned int r,aux;
   for(i = n-1; i >= 1; i--)
   {
       r = R[n-i] % (i+1);
       aux = perm[i];
       perm[i] = perm[r];
       perm[r] = aux;
   }
   for(i = 0; i<n; i++)
    Lprim[perm[i]] = L[i];

   ///SUBSTITUTIE
   unsigned char sv_octet[3], r_octet[3];
   sv_octet[0] =  SV        & 255;
   sv_octet[1] = (SV >> 8)  & 255;
   sv_octet[2] = (SV >> 16) & 255;

   r_octet[0] =  R[n]        & 255;
   r_octet[1] = (R[n] >> 8)  & 255;
   r_octet[2] = (R[n] >> 16) & 255;

   Lprim[0].R = sv_octet[2] ^ Lprim[0].R ^ r_octet[2];
   Lprim[0].G = sv_octet[1] ^ Lprim[0].G ^ r_octet[1];
   Lprim[0].B = sv_octet[0] ^ Lprim[0].B ^ r_octet[0];

   for(i=1;i<n;i++)
   {
       r_octet[0] =  R[n+i]        & 255;
       r_octet[1] = (R[n+i] >> 8)  & 255;
       r_octet[2] = (R[n+i] >> 16) & 255;

       Lprim[i].R = Lprim[i-1].R ^ Lprim[i].R ^ r_octet[2];
       Lprim[i].G = Lprim[i-1].G ^ Lprim[i].G ^ r_octet[1];
       Lprim[i].B = Lprim[i-1].B ^ Lprim[i].B ^ r_octet[0];
   }

   salvare_imagine (nume_fisier_sursa, nume_fisier_destinatie, Lprim);

   free(R); free(L); free(Lprim); free(perm);
}

void decriptare(char *nume_fisier_sursa, char *nume_fisier_destinatie, char *cheie_secreta)
{
   ///CHEIA
   unsigned int R0, SV;
   FILE *key;
   key = fopen(cheie_secreta, "r");
   if(key == NULL)
   {
       printf("Fisierul unde se afla cheia nu a putut fi deschis");
       return;
   }
   fscanf(key, "%u", &R0); fscanf(key, "%u", &SV);
   fclose(key);


   FILE *fin;
   fin = fopen(nume_fisier_sursa, "rb");
   unsigned int latime_img, inaltime_img, n;

   if(fin == NULL)
   {
       printf("Nu am putut gasi imaginea");
       return;
   }

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img,   sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);
   fclose(fin);

   n = latime_img * inaltime_img;

   pixel *L, *Lprim;
   unsigned int *R, i;

   R = XORSHIFT32 (R0,n);
   L = incarca_imagine (nume_fisier_sursa);

   ///SUBSTITUTIE
   unsigned char sv_octet[3], r_octet[3];

   for(i = n-1; i >= 1; i--)
   {
       r_octet[0] =  R[n+i]        & 255;
       r_octet[1] = (R[n+i] >> 8)  & 255;
       r_octet[2] = (R[n+i] >> 16) & 255;

       L[i].R = L[i-1].R ^ L[i].R ^ r_octet[2];
       L[i].G = L[i-1].G ^ L[i].G ^ r_octet[1];
       L[i].B = L[i-1].B ^ L[i].B ^ r_octet[0];
   }

   sv_octet[0] =  SV        & 255;
   sv_octet[1] = (SV >> 8)  & 255;
   sv_octet[2] = (SV >> 16) & 255;

   r_octet[0] =  R[n]        & 255;
   r_octet[1] = (R[n] >> 8)  & 255;
   r_octet[2] = (R[n] >> 16) & 255;

   L[0].R = sv_octet[2] ^ L[0].R ^ r_octet[2];
   L[0].G = sv_octet[1] ^ L[0].G ^ r_octet[1];
   L[0].B = sv_octet[0] ^ L[0].B ^ r_octet[0];

   Lprim=(pixel *) calloc (n, 3);
   unsigned int *perm     = calloc (n, sizeof(unsigned int));
   unsigned int *inv_perm = calloc (n, sizeof(unsigned int));

   ///PERMUTARE
   for(i = 0; i < n; i++) perm[i] = i;
   unsigned int r,aux;
   for(i = n-1; i >= 1; i--)
   {
       r = R[n-i] % (i+1);
       aux = perm[i];
       perm[i] = perm[r];
       perm[r] = aux;
   }
   for(i = 0; i < n; i++) inv_perm[perm[i]] = i;
   for(i = 0; i < n; i++) Lprim[inv_perm[i]] = L[i];

   salvare_imagine (nume_fisier_sursa, nume_fisier_destinatie, Lprim);

   free(R); free(L); free(Lprim); free(perm); free(inv_perm);
}

void chi_patrat(char *nume_fisier)
{
    FILE *fin;
    fin = fopen(nume_fisier, "rb");
    if(fin == NULL)
    {
        printf("Nu s-a putut deschide imaginea in functia chi_patrat");
        return;
    }

    unsigned int latime_img, inaltime_img;

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

    unsigned int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    fseek(fin, 54, SEEK_SET);

    unsigned int i,j;
    unsigned char red, green, blue;

    unsigned int *fr_R = calloc (256, sizeof(unsigned int));
    unsigned int *fr_G = calloc (256, sizeof(unsigned int));
    unsigned int *fr_B = calloc (256, sizeof(unsigned int));

	for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++){
            fread(&blue,  1, 1, fin); fr_B[blue]++;
            fread(&green, 1, 1, fin); fr_G[green]++;
            fread(&red,   1, 1, fin); fr_R[red]++;}
        fseek(fin, padding, SEEK_CUR);
    }

    double frecventa_teoretica = (inaltime_img*latime_img)/256.0;
    double chi_blue=0, chi_green=0, chi_red=0;

    for(i=0;i<256;i++)
    {
        chi_blue  =  chi_blue  + ((fr_B[i]-frecventa_teoretica)*(fr_B[i]-frecventa_teoretica))/frecventa_teoretica;
        chi_green =  chi_green + ((fr_G[i]-frecventa_teoretica)*(fr_G[i]-frecventa_teoretica))/frecventa_teoretica;
        chi_red   =  chi_red   + ((fr_R[i]-frecventa_teoretica)*(fr_R[i]-frecventa_teoretica))/frecventa_teoretica;
    }

    printf("%0.2f %0.2f %0.2f\n", chi_red, chi_green, chi_blue);

    fclose(fin);
}

void grayscale_image(char* nume_fisier_sursa, char* nume_fisier_destinatie)
{
   FILE *fin, *fout;
   unsigned int latime_img, inaltime_img;
   unsigned char pRGB[3], aux;

   fin = fopen(nume_fisier_sursa, "rb");
   if(fin == NULL)
   	{
   		printf("Nu am gasit imaginea sursa din care citesc, in functia grayscale");
   		return;
    }

    fout = fopen(nume_fisier_destinatie, "wb+");

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);

	fseek(fin,0,SEEK_SET);
	unsigned char c;
	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);

	int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

	fseek(fout, 54, SEEK_SET);
	int i,j;
	for(i = 0; i < inaltime_img; i++)
	{
		for(j = 0; j < latime_img; j++)
		{
			fread(pRGB, 3, 1, fout);
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
        	fseek(fout, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, fout);
        	fflush(fout);
		}
		fseek(fout,padding,SEEK_CUR);
	}
	fclose(fout);
	printf("Imaginea a fost facuta gri\n");
}

void contur (pixel *I, pixel C, fereastra f, unsigned int inaltime, unsigned int latime, unsigned int latimeI)
{
    unsigned int i, j;
    for(i = 0; i < inaltime; i++)
    {
        I[(f.i + i) * latimeI + f.j].B = C.B;
        I[(f.i + i) * latimeI + f.j].G = C.G;
        I[(f.i + i) * latimeI + f.j].R = C.R;

        I[(f.i + i) * latimeI + f.j + latime - 1].B = C.B;
        I[(f.i + i) * latimeI + f.j + latime - 1].G = C.G;
        I[(f.i + i) * latimeI + f.j + latime - 1].R = C.R;
    }


    for(j = 0; j < latime; j++)
    {
        I[f.i * latimeI + f.j + j].B = C.B;
        I[f.i * latimeI + f.j + j].G = C.G;
        I[f.i * latimeI + f.j + j].R = C.R;

        I[(f.i + inaltime - 1) * latimeI + f.j + j].B = C.B;
        I[(f.i + inaltime - 1) * latimeI + f.j + j].G = C.G;
        I[(f.i + inaltime - 1) * latimeI + f.j + j].R = C.R;
    }
}

double valoare_intensitate_sablon (char *nume_sablon, unsigned int inaltime_img, unsigned int latime_img)
{
    FILE *fin;
    fin = fopen (nume_sablon, "rb");
    if(fin == NULL)
    {
        printf("Eroare la deschiderea sablonului");
        return;
    }

    unsigned int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    fseek(fin, 54, SEEK_SET);
    unsigned int i, j;
    double S_barat = 0;
    unsigned char c;

	for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            fread(&c, 1, 1, fin);
            fseek(fin, 2, SEEK_CUR);
            S_barat += c;
        }
        fseek(fin, padding, SEEK_CUR);
    }
    fclose(fin);

    S_barat = S_barat/(double)(latime_img * inaltime_img);

    return S_barat;
}

double valoare_intensitate_fereastra (pixel *I, fereastra f, unsigned int inaltime_img, unsigned int latime_img, unsigned int latimeI)
{
    unsigned int i, j;
    double F_barat = 0;

    for(i = 0; i < inaltime_img; i++)
        for(j = 0; j < latime_img; j++)
            F_barat += I[(f.i + i) * latimeI + f.j + j].R;

    F_barat = F_barat/(double)(latime_img * inaltime_img);

    return F_barat;
}

double deviatie_sablon (char *nume_sablon, unsigned int inaltime_img, unsigned int latime_img)
{
    double S_barat = valoare_intensitate_sablon (nume_sablon, inaltime_img, latime_img);

    FILE *fin;
    fin = fopen (nume_sablon, "rb");
    if(fin == NULL)
    {
        printf("Eroare la deschiderea sablonului");
        return;
    }

    unsigned int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    fseek(fin, 54, SEEK_SET);
    unsigned int i, j;
    double S_deviatie = 0;
    unsigned char c;

	for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            fread(&c, 1, 1, fin);
            fseek(fin, 2, SEEK_CUR);
            S_deviatie += (c - S_barat) * (c - S_barat);
        }
        fseek(fin, padding, SEEK_CUR);
    }
    fclose(fin);

    S_deviatie = S_deviatie/(double)(latime_img * inaltime_img - 1);
    S_deviatie = sqrt (S_deviatie);

    return S_deviatie;
}

double deviatie_fereastra (pixel *I, fereastra f, unsigned int inaltime_img, unsigned int latime_img, unsigned int latimeI)
{

    unsigned int i, j;
    double F_barat = valoare_intensitate_fereastra (I, f, inaltime_img, latime_img, latimeI);
    double F_deviatie = 0;

    for(i = 0; i < inaltime_img; i++)
        for(j = 0; j < latime_img; j++)
            F_deviatie += (I[(f.i + i) * latimeI + f.j + j].R - F_barat) * (I[(f.i + i) * latimeI + f.j + j].R - F_barat);

    F_deviatie = F_deviatie/(double)(latime_img * inaltime_img - 1);
    F_deviatie = sqrt (F_deviatie);

    return F_deviatie;
}

double corelatie (char *nume_sablon, pixel *I, fereastra f, unsigned int inaltime_img, unsigned int latime_img, unsigned int latimeI, double S_barat, double S_deviatie)
{
    double F_barat, F_deviatie, corelatiee = 0 ;
    F_barat = valoare_intensitate_fereastra (I, f, inaltime_img, latime_img, latimeI);
    F_deviatie = deviatie_fereastra (I, f, inaltime_img, latime_img, latimeI);

    FILE *fin;
    fin = fopen (nume_sablon, "rb");
    if(fin == NULL)
    {
        printf("Eroare la deschiderea sablonului");
        return;
    }

    unsigned int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    fseek(fin, 54, SEEK_SET);
    unsigned int i, j;
    unsigned char c;

	for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            fread(&c, 1, 1, fin);
            fseek(fin, 2, SEEK_CUR);
            corelatiee += ((I[(f.i + i) * latimeI + f.j + j].R - F_barat) * (c - S_barat));
        }
        fseek(fin, padding, SEEK_CUR);
    }
    fclose(fin);

    corelatiee = corelatiee / (latime_img * inaltime_img * F_deviatie * S_deviatie);

    return corelatiee;
}

unsigned int suprapunere (fereastra f1, fereastra f2, unsigned int inaltime, unsigned int latime)
{
    unsigned int i = abs(f1.i - f2.i);
    unsigned int j = abs(f1.j - f2.j);

    if(i >= inaltime && j >= latime) return 1;

    unsigned int arie_intersectie;
    double supra;
    arie_intersectie = (inaltime - i) * (latime - j);
    supra = (double)arie_intersectie / (double)(2*inaltime*latime - arie_intersectie);

    if(supra <= 0.2) return 1;

  return 0;
}

int cmp (const void *d_1, const void *d_2)
{
    detectie *x = (detectie*)d_1;
    detectie *y = (detectie*)d_2;

    if(x->corelatiee < y->corelatiee)
        return 1;
    return -1;
}

pixel* incarca_imagine2 (char *nume_fisier_sursa)
{
   FILE *fin;
   unsigned int latime_img, inaltime_img;
   pixel *L;

   fin = fopen(nume_fisier_sursa, "rb");
   if(fin == NULL)
   {
       printf("Nu am putut gasii imaginea");
       return;
   }

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img,   sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);


   unsigned int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;



    fseek(fin, 54, SEEK_SET);
    unsigned int i,j;

    L = (pixel *) calloc (latime_img*inaltime_img, 3);
    if(L == NULL)
    {
        printf("Eroare la alocare");
        return;
    }

	for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++){
            fread(&(L[i * latime_img + j].B), 1, 1, fin);
            fread(&(L[i * latime_img + j].G), 1, 1, fin);
            fread(&(L[i * latime_img + j].R), 1, 1, fin);}
        fseek(fin, padding, SEEK_CUR);
    }

    fclose(fin);
    return L;
}

void salvare_imagine2 (char *nume_fisier_sursa, pixel *C)
{
    FILE *fin;
    unsigned int latime_img, inaltime_img, i, j;

    fin = fopen(nume_fisier_sursa, "rb+");

    if(fin == NULL)
    {
        printf("Nu am putut gasii imaginea");
        return 0;
    }

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);

   unsigned int padding;
   unsigned char padd = 0x00;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

    fseek(fin,54,SEEK_SET);


	for(i = 0; i < inaltime_img; i++)
    {
        for(j = 0; j < latime_img; j++)
        {
            fwrite(&(C[i * latime_img + j].B), 1, 1, fin);
            fwrite(&(C[i * latime_img + j].G), 1, 1, fin);
            fwrite(&(C[i * latime_img + j].R), 1, 1, fin);
            fflush(fin);
        }
        for(int pd = 0; pd < padding; pd++){
        fwrite(&padd, 1, 1, fin); fflush(fin);}
    }

    fclose(fin);
    printf("Imaginea decriptata a fost modificata\n\n");
}

void template_maching (char *nume_imagine, char *nume_sablon, double ps, pixel color, unsigned int *nr, detectie *D)
{
    pixel *I;
    I = incarca_imagine2 (nume_imagine);

    unsigned int latimeI, latime_img, inaltime_img, inaltimeI;
    FILE *ffff ,*fin;
    ffff = fopen (nume_imagine, "rb");
    if(ffff == NULL)
    {
        printf("Eroare la deschidere fisier in tepmlate_maching");
        return;
    }
    fseek(ffff, 18, SEEK_SET);
    fread(&latimeI,   sizeof(unsigned int), 1, ffff);
    fread(&inaltimeI, sizeof(unsigned int), 1, ffff);
    fclose(ffff);

    fin = fopen(nume_sablon, "rb");
    if(fin == NULL)
    {
       printf("Nu am putut gasii imaginea");
       return;
    }

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img,   sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);
    fclose(fin);


    double S_barat, S_deviatie, corelation;
    S_barat    = valoare_intensitate_sablon (nume_sablon, inaltime_img, latime_img);
    S_deviatie = deviatie_sablon (nume_sablon, inaltime_img, latime_img);

    unsigned int a, b; fereastra ff;
    for(int a = 0; a < (inaltimeI - inaltime_img + 1); a++)
    {
        ff.i = a;
        for(int b = 0; b < (latimeI - latime_img + 1); b++)
        {
            ff.j = b;
            corelation = corelatie (nume_sablon, I, ff, inaltime_img, latime_img, latimeI, S_barat, S_deviatie);
            if(corelation >= ps)
            {
                D[(*nr)].corelatiee = corelation;
                D[(*nr)].culoare = color;
                D[(*nr)].exista = 1;
                D[(*nr)].fe.i = ff.i;
                D[(*nr)].fe.j = ff.j;
                (*nr)++;
            }
        }
    }


    free(I);
}

void eliminare_non_maxime (detectie *D, unsigned int numar_detectii, char *nume_sablon)
{
    qsort (D, numar_detectii, sizeof(detectie), cmp);

    unsigned int latime_img, inaltime_img;
    FILE *fin;
    fin = fopen(nume_sablon, "rb");
    if(fin == NULL)
    {
       printf("Nu am putut gasii imaginea");
       return;
    }

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img,   sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);
    fclose(fin);

    unsigned int a, b;
    for(a = 0; a < numar_detectii - 1; a++)
        for(b = a + 1; b < numar_detectii; b++)
            if(D[a].exista == 1 && D[b].exista == 1)
               if((D[a].fe.i - D[b].fe.i < inaltime_img-2 || D[a].fe.i - D[b].fe.i > -inaltime_img+2) && (D[a].fe.j - D[b].fe.j < latime_img-2 || D[a].fe.j - D[b].fe.j > -latime_img+2))
                D[b].exista = 0;
}


int main()
{

    char nume_imagine_sursa[101], nume_imagine_destinatie[101], nume_imagine_decriptata[101], nume_fisier_cheie[101];
    char nume_sabloane[] = "nume_sabloane.txt", nume_imagine_gray[101], culori[] = "Culori.txt";

    printf("Numele fisierului care contine imaginea pe care o vom cripta: ");
    fgets(nume_imagine_sursa, 101, stdin);
    nume_imagine_sursa[strlen(nume_imagine_sursa) - 1] = '\0';

    printf("Numele fisierului care va contine imaginea criptata: ");
    fgets(nume_imagine_destinatie, 101, stdin);
    nume_imagine_destinatie[strlen(nume_imagine_destinatie) - 1] = '\0';

    printf("Numele fisierului care contine cheia: ");
    fgets(nume_fisier_cheie, 101, stdin);
    nume_fisier_cheie[strlen(nume_fisier_cheie) - 1] = '\0';

    /// CRIPTARE !!!
    criptare(nume_imagine_sursa, nume_imagine_destinatie, nume_fisier_cheie);

    printf("Numele fisierului care va contine imaginea decriptata: ");
    fgets(nume_imagine_decriptata, 101, stdin);
    nume_imagine_decriptata[strlen(nume_imagine_decriptata) - 1] = '\0';

    ///DECRIPTARE !!!
    decriptare(nume_imagine_destinatie, nume_imagine_decriptata, nume_fisier_cheie);

    printf("Tetul chi patrat:\n");
    ///TESTUL CHI !!!
    chi_patrat(nume_imagine_sursa);
    chi_patrat(nume_imagine_destinatie);

    printf("Numele fisierului care va contine imaginea gri: ");
    fgets(nume_imagine_gray, 101, stdin);
    nume_imagine_gray[strlen(nume_imagine_gray) - 1] = '\0';
    ///CREEAM IMAGINILE GRI
    grayscale_image(nume_imagine_decriptata, nume_imagine_gray);





    FILE *cul, *sab;
    cul = fopen (culori, "r");
    if(cul == NULL)
    {
        printf("Fisierul in care se afla culorile nu a putut fi deschis");
        return;
    }

    sab = fopen (nume_sabloane, "r");
    if(sab == NULL)
    {
        printf("Fisierul in care se afla numele sabloanelor nu a putut fi deschis");
        return;
    }




    pixel C;
    char nume_sablon[40];
    unsigned int numar_detectii = 0, i;
    detectie *D = calloc (3500, sizeof(detectie));

    if(D == NULL)
    {
        printf("Eroare la alocare");
        return;
    }

    for(i = 0; i < 10; i++)
    {
        fscanf(sab,"%s", nume_sablon);
        fscanf(cul,"%u", &(C.R));
        fscanf(cul,"%u", &(C.G));
        fscanf(cul,"%u", &(C.B));
        template_maching (nume_imagine_gray, nume_sablon, 0.5, C, &numar_detectii, D);
    }

    D = (detectie *) realloc (D, numar_detectii * sizeof(detectie));
    if(D == NULL)
    {
        printf("Eroare la alocare");
        return;
    }
    fclose(sab); fclose(cul);

    eliminare_non_maxime (D, numar_detectii, nume_sablon);


    pixel *I_final;
    I_final = incarca_imagine2 (nume_imagine_sursa);


    for(i = 0; i < numar_detectii; i++)
        if(D[i].exista == 1)
            contur (I_final, D[i].culoare, D[i].fe, 15, 11, 500);


     salvare_imagine2 (nume_imagine_decriptata, I_final);



     free(D); free(I_final);

     return 0;
}
