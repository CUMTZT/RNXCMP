/***************************************************************************
 * ��飺 RNX ��ʽѹ���㷨ʵ��
 * �÷���
          RNX2CRX [file] [-] [-f] [-e # of epochs] [-s] [-d] [-h]
            stdin and stdout are used if input file name is not given.
            -       : output to stdout
            -f      : force overwrite of output file
            -e #    : initialize the compression operation at every # epochs
                        When some part of the Compact RINEX file is lost, the data
                        can not be recovered thereafter until all the data arc are
                        initialized for differential operation. This option may be used to
                        increase chances to recover parts of data by using an option of
                        CRX2RNX(ver. 4.0 or after) with cost of increase of file size.
            -s      : warn and skip strange epochs (default: stop with error status)
            -d      : delete the input file if conversion finishes without errors
                      (i.e. exit code = 0 or 2).
                      This option does nothing if stdin is used for the input.
            -h      : display help message
****************************************************************************/

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

/***************************************************************************
* ��飺����汾��
* ��;���������,���㷨�޹أ�����
****************************************************************************/

#define VERSION  "ver.4.1.0"

/***************************************************************************
* ��飺����˳���
* ��;������ʶ������˳�״̬�����㷨�޹أ�����
****************************************************************************/
#define EXIT_WARNING 2

#ifndef EXIT_FAILURE
#define EXIT_FAILURE 1
#endif

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#endif

/***************************************************************************
* ��飺�ڴ������
* ��;������ÿ��ѭ��������һ��ѭ��������
*       FLUSH_BUFF �ȴ�ӡ��������
*       CLEAR_BUFF ֱ������
****************************************************************************/
#define FLUSH_BUFF printf("%s",top_buff), *(p_buff = top_buff)='\0'
#define CLEAR_BUFF *(p_buff = top_buff) = '\0'

/***************************************************************************
* ��飺һЩ�����궨��
* ��;������һЩѹ���㷨�汾�š��ڴ�������󳤶�
*       �汾����Ҫ�������
*       �ڴ�������󳤶���Ҫ���������ڴ�ʱȷ�����ֵ
****************************************************************************/
#define CRX_VERSION1 "1.0"    /* CRINEX version for RINEX 2.x */
#define CRX_VERSION2 "3.0"    /* CRINEX version for RINEX 3.x */
#define PROGNAME "RNX2CRX"
#define MAXSAT    100         /* Maximum number of satellites observed at one epoch */
#define MAXTYPE   100         /* Maximum number of data types for a GNSS system */
#define MAXCLM   2048         /* Maximum columns in one line   (>MAXTYPE*19+3)  */
#define MAX_BUFF_SIZE 204800  /* Maximum size of output buffer (>MAXSAT*(MAXTYPE*19+4)+60 */
#define ARC_ORDER 3           /* order of difference to take    */

/* define data structure for fields of clock offset and observation records */
/* Those data will be handled as integers after eliminating decimal points.  */
/* Since their size may exceeds range between LONG_MIN and LONG_MAX,         */
/* they will be read with being divided properly into upper and lower digits */

/***************************************************************************
* ��飺ʱ���ʽ�ṹ��
* ��;�����ڴ洢ʱ�䣬�������ݿ��ܳ���long�����ݷ�Χ����˽�ʱ�������з�
****************************************************************************/
typedef struct clock_format {
    long u[ARC_ORDER + 1];      //��λ
    long l[ARC_ORDER + 1];      //��λ
} clock_format;

/***************************************************************************
* ��飺���ݸ�ʽ�ṹ��
* ��;�����ڴ洢���ݣ��������ݿ��ܳ���long�����ݷ�Χ����˽�ʱ�������з�
****************************************************************************/
typedef struct data_format {
    long u[ARC_ORDER + 1];      //��λ
    long l[ARC_ORDER + 1];      //��λ
    int  order;
} data_format;

/***************************************************************************
* ��飺ȫ�ֱ�������
* ��;������ȫ�ֱ����Ķ��弰��ʼ�����ⲿ�ֶ���ı���һ��Ҫ���ÿһ������
*       �ĺ��壬�����ڶ������㷨�����
****************************************************************************/
long ep_count = 0;
long ep_reset = 0;
long nl_count = 0;
int rinex_version;          /* =2, 3 or 4 */
int nsat, ntype, ntype_gnss[UCHAR_MAX], ntype_record[MAXSAT], clk_order = -1;
int exit_status = EXIT_SUCCESS;
int skip_strange_epoch = 0; /* default : stop with error */
int delete_if_no_error = 0; /* default : not delete */
int n_infile = 0;           /* number of input file (must be 0 or 1) */


clock_format clk1 = { {0,0,0,0},{0,0,0,0} }; //ʱ��ṹ��������ڴ洢ʱ��
clock_format clk0 = { {0,0,0,0},{0,0,0,0} }; //ʱ��ṹ��������ڴ洢ʱ��
data_format dy0[MAXSAT][MAXTYPE], dy1[MAXSAT][MAXTYPE]; //���ݽṹ�������������ʱ��
char flag0[MAXSAT][MAXTYPE * 2], flag[MAXSAT][MAXTYPE * 2];
char out_buff[MAX_BUFF_SIZE] = { 'x','\0' }; //���Buffer�洢���飬��λΪ��x��������λ�ÿ�
char* top_buff = &out_buff[1], * p_buff;        //���������λ��ȻΪ��x����ʣ�µ�λ����Ϊ���ɵ�����
char infile[MAXCLM];                        //�����ļ����洢����
constexpr size_t C1 = sizeof("");               //һ���ַ����ȣ����߿��ܿ��ǵ�������ַ����ȵ�Ӱ�죬��������ַ������ж���
constexpr size_t C2 = sizeof(" ");              //�����ַ�����
constexpr size_t C3 = sizeof("  ");             //�����ַ�����
constexpr size_t C14 = sizeof("             "); //ʮ�ĸ��ַ�����

char oldline[MAXCLM] = { '&','\0' };
int nsat_old = 0;

/***************************************************************************
* ��飺����������
* ��;�������㷨һЩ����������������main�������ã����嶨���ں��
****************************************************************************/
void parse_args(int argc, char* argv[]);//��������Ĳ��������ݲ�����ʼ��һЩȫ�ֱ�����Դ�ļ����ƣ�ѹ���ɹ����Ƿ�ɾ��Դ�ļ��Ȳ������ã�
void header(void);//���������ļ���ͷ������ͷ��ʼ��һЩȫ�ֱ�����Դ�ļ���RNX�ļ��İ汾�ȱ�����
int  get_next_epoch(char* p_line);
void skip_to_next(char* p_line);
void initialize_all(char* oldline, int* nsat_old, int count);
void put_event_data(char* p_line);
void read_clock(char* line, int shift_cl);
void process_clock(void);
int  set_sat_table(char* p_new, char* p_old, int nsat_old, int* sattbl);
int  read_more_sat(int n, char* p);
void data(int* sattbl);
char* strdiff(const char* s1, char* s2, char* ds);
int  ggetline(data_format* py1, char* flag, char* sat_id, int* ntype_rec);
void read_value(char* p, long* pu, long* pl);
void take_diff(data_format* py1, data_format* py0);
void putdiff(long dddu, long dddl);
void put_clock(long du, long dl, int clk_order);
int  read_chk_line(char* line);
void error_exit(int error_no, char* string);
void no_error_exit();


/***************************************************************************
* ��飺������
* ��;���������������㷨�������ڣ�Ҳ����˵ʹ�������������ó���ʱ��argc��Ϊ
*       ����ʱ�Ĳ���������argv�Ǹ��ַ������飬�������argc���ַ������������
*       ��������Ϊ��.\RNX2CRX.exe test.rnx��,��argc��ֵΪ2��argv[0]Ϊ
*        ��.\RNX2CRX.exe����argv[1]Ϊtest.rnx����������Դ����ơ�
****************************************************************************/
int main(int argc, char* argv[]) {

    char newline[MAXCLM];
    char dummy[2] = { '\0','\0' };
    char* p, * p_event, * p_nsat, * p_satlst, * p_satold, * p_clock;
    //��¼�������ϸ�epoch line��λ��
    int sattbl[MAXSAT];

    int i, j, shift_clk;

    parse_args(argc, argv);//��������,���㷨�޹أ��ɲ���ע

    //����GNSS�ĳ�ֵΪ-1��������ļ���û�ж��壬��ʹ��-1
    for (i = 0; i < UCHAR_MAX; i++) ntype_gnss[i] = -1;  /** -1 unless GNSS type is defined **/

    //�����ļ�ͷ
    header();

    //�����볤��ΪMAXCLM�����飬Ȼ����ݰ汾�Ž�ָ��ָ��Ӧ���������ݵ�λ��
    if (rinex_version == 2) {
        p_event = &newline[28];  //ָ��event flag
        p_nsat = &newline[29];  //ָ��n_sat
        p_satlst = &newline[32];  //ָ��satellite list
        p_satold = &oldline[32];  //ָ����һ�е�n_sat
        p_clock = &newline[68];  //ָ��clock offset data
        shift_clk = 1;
    }
    else {
        p_event = &newline[31];
        p_nsat = &newline[32];
        p_satlst = &newline[41];
        p_satold = &oldline[41];
        p_clock = &newline[41];
        shift_clk = 4;
    }

/***************************************************************************
* ��飺��������ѭ������
* ע�ͣ�forѭ��
*       ����ִ��CLEAR_BUFF,Ȼ��ִ��forѭ�����ڵĴ��룬ִ�����ִ��FLUSH_BUFF
*       ��forѭ��ֻ����exitʱ�˳�
****************************************************************************/
    for (CLEAR_BUFF;; FLUSH_BUFF) {
        //���ص�ʱ0���ļ������������˳�����
        if (!get_next_epoch(newline)) no_error_exit();

        //���������p_event��Ӧ��ֵ����1������Ҫ���event data
        if (atoi(strncpy(dummy, p_event, C1)) > 1) {
            put_event_data(newline);
            initialize_all(oldline, &nsat_old, 0);
            continue;
        }


        //��������ˡ�-e������ֵ����0�����ep_count �����˸�ֵ�������ò���
        if (ep_reset > 0 && ++ep_count > ep_reset) initialize_all(oldline, &nsat_old, 1);


        //�жϸ�����û��ʱ�����ݣ��еĻ���ȡ�������˼��������ļ���û����������ݵ������ļ���
        if (strchr(newline, '\0') > p_clock) {
            read_clock(p_clock, shift_clk);        /**** read clock offset ****/
        }
        else {
            clk_order = -1;                       /*** reset data arc for clock offset ***/
        }

        //��ȡ�м�������
        nsat = atoi(p_nsat);
        if (nsat > MAXSAT) error_exit(8, newline);
        if (nsat > 12 && rinex_version == 2) read_more_sat(nsat, p_satlst);  /*** read continuation lines ***/

        //�����������������ж�ȡ���ݣ����������浽buffer���C01C04C08C11C12C19C21C22��
        for (i = 0, p = p_satlst; i < nsat; i++, p += 3) {
            if (ggetline(dy1[i], flag[i], p, &ntype_record[i])) {
                CLEAR_BUFF;
                exit_status = EXIT_WARNING;
                continue;
            }
        }
        *p = '\0';   

        //������һ�е�����˳��ͱ�������˳�����ɱ�����������һ�еĶ�Ӧλ��
        if (set_sat_table(p_satlst, p_satold, nsat_old, sattbl)) {
            CLEAR_BUFF;
            exit_status = EXIT_WARNING;
            continue;
        }


        //��ʱ���е���ͷ�Ѿ�������ϣ�ͨ����������������һ���������Աȣ����ɱ���Ӧ�������ѹ������
        p_buff = strdiff(oldline, newline, p_buff);

        //û�ж�Ӧ���ݣ���ʱû��
        if (clk_order > -1) {
            if (clk_order > 0) process_clock();            /**** process clock offset ****/
            put_clock(clk1.u[clk_order], clk1.l[clk_order], clk_order);
        }
        else {
            *p_buff++ = '\n';
        }

        data(sattbl); *p_buff = '\0';
        /**************************************/
        /**** save current epoch to buffer ****/
        /**************************************/
        nsat_old = nsat;
        sprintf(oldline, "%s", newline);
        clk0 = clk1;
        for (i = 0; i < nsat; i++) {
            strcpy(flag0[i], flag[i]);
            for (j = 0; j < ntype_record[i]; j++) dy0[i][j] = dy1[i][j];
        }
    }
}
/*---------------------------------------------------------------------*/
void parse_args(int argc, char* argv[]) {
    char* p;
    char outfile[MAXCLM];//����ļ�������XX.crx��
    char* progname;//�����
    int force = 0, help = 0;
    int nfout = 0;  /*** =0 default output file name ***/
    /*** =1 standard output          ***/
    FILE* ifp;

    progname = argv[0];//�������Ϊ��һ������

    //����һ���������������ȥ��
    argc--; argv++;
    //ѭ���������еĲ���
    for (; argc > 0; argc--, argv++) {
        //������һ���ַ����ǡ�-�����Ǳ�Ȼ���ļ���
        if ((*argv)[0] != '-') {
            //�������ļ��������infile���������
            strncpy(infile, *argv, C1 * MAXCLM);
            n_infile++;
        }
        //������-��
        else if (strcmp(*argv, "-") == 0) {
            nfout = 1;                 /* output to standard output */
        }
        //������-f��,�������ļ�����������
        else if (strcmp(*argv, "-f") == 0) {
            force = 1;                 /* overwrite if the output file exists */
        }
        //������-d��,���ѹ��û�д���ɾ��Դ�ļ�
        else if (strcmp(*argv, "-d") == 0) {
            delete_if_no_error = 1;    /* delete the original file if
                                          no error in the conversion */
        }
        //������-s��,������������
        else if (strcmp(*argv, "-s") == 0) {
            skip_strange_epoch = 1;
        }
        else if (strcmp(*argv, "-e") == 0) {
            argc--; argv++;
            sscanf(*argv, "%ld", &ep_reset);
        }
        else if (strcmp(*argv, "-h") == 0) {
            help = 1;
        }
        else {
            help = 1;
        }
    }

    if (strlen(infile) == MAXCLM) error_exit(14, infile);
    if (help == 1 || n_infile > 1 || n_infile < 0) error_exit(1, progname);
    if (n_infile == 0) return;  /*** stdin & stdout will be used if input file name is not given ***/

    /***********************/
    /*** open input file ***/
    /***********************/
    p = strrchr(infile, '.');
    if (p == NULL || *(p + 4) != '\0'
        || (toupper(*(p + 3)) != 'O'
            && strcmp(p + 1, "RNX") != 0
            && strcmp(p + 1, "rnx") != 0)
        ) error_exit(4, p);

    if ((ifp = fopen(infile, "r")) == NULL) error_exit(5, infile);

    /************************/
    /*** open output file ***/
    /************************/
    if (nfout == 0) {
        strcpy(outfile, infile);
        p = strrchr(outfile, '.');
        if (*(p + 3) == 'o') { *(p + 3) = 'd'; }
        else if (*(p + 3) == 'O') { *(p + 3) = 'D'; }
        else if (strcmp(p + 1, "rnx") == 0) { strcpy((p + 1), "crx"); }
        else if (strcmp(p + 1, "RNX") == 0) { strcpy((p + 1), "CRX"); }

        if ((freopen(outfile, "r", stdout)) != NULL && force == 0) {
            fprintf(stderr, "The file %s already exists. Overwrite?(n)", outfile);
            if (getchar() != 'y') exit(EXIT_SUCCESS);
        }
        freopen(outfile, "w", stdout);
    }
    fclose(ifp);
    freopen(infile, "r", stdin);
}
/*---------------------------------------------------------------------*/
void header(void) {
    //�������ݴ洢����
    char line[MAXCLM], line2[41], timestring[20];
    //����ʱ��ṹ�����
    time_t tc = time(NULL);
    struct tm* tp;

    //��ȡ��ǰ��UTCʱ�䣬����ǰ����ʱ��ȥ8Сʱ��ʾ�����ݣ���02-Apr-24 02:01��
    if ((tp = gmtime(&tc)) == NULL) tp = localtime(&tc);
    strftime(timestring, C1 * 20, "%d-%b-%y %H:%M", tp);

    //���rnx�ļ��ĵ�һ���Ƿ������ȷ
    read_chk_line(line);
    if (strncmp(&line[60], "RINEX VERSION / TYPE", C1 * 20) != 0 ||
        strncmp(&line[20], "O", C1) != 0) error_exit(15, line);

    //���ݰ汾�����crx�ļ��ĵ�һ��
    rinex_version = atoi(line);
    if (rinex_version == 2) { printf("%-20.20s", CRX_VERSION1); }
    else if (rinex_version == 3 || rinex_version == 4) { printf("%-20.20s", CRX_VERSION2); }
    else { error_exit(15, line); }
    printf("%-40.40s%-20.20s\n", "COMPACT RINEX FORMAT", "CRINEX VERS   / TYPE");

    //���crx�ļ��ڶ���
    sprintf(line2, "%s %s", PROGNAME, VERSION);
    printf("%-40.40s%-20.20sCRINEX PROG / DATE\n", line2, timestring);
    //��rnx�ĵ�һ��ֱ�������crx�ĵ�����
    printf("%s\n", line);

    //ѭ����ȡ��ֱ��������END OF HEADER��
    do {
        read_chk_line(line);

        //���������������crx�ļ���
        printf("%s\n", line);

        //����������
        if (strncmp(&line[60], "# / TYPES OF OBSERV", C1 * 19) == 0 && line[5] != ' ') {
            ntype = atoi(line);                                        /** for RINEX2 **/
        }
        else if (strncmp(&line[60], "SYS / # / OBS TYPES", C1 * 19) == 0) { /** for RINEX3 **/
            if (line[0] != ' ') ntype_gnss[(unsigned int)line[0]] = atoi(&line[3]);
            if (ntype_gnss[(unsigned int)line[0]] > MAXTYPE) error_exit(16, line);
        }
    } while (strncmp(&line[60], "END OF HEADER", C1 * 13) != 0);
}

/***************************************************************************
* ��飺��ȡ��һ����Ч��
* ע�ͣ�������и�ʽ����ȷ��������ֱ����ȡ����ȷ��һ��
* ����ֵ��
*       0����ȡ�����ļ����һ�У���û�л�ȡ����Ӧ�����ݣ������㷨
*       1����ȡ������ȷһ�У��㷨����
*       2����ȡ���˴����һ��
****************************************************************************/
int  get_next_epoch(char* p_line) {
    char* p;

    nl_count++;
    //��ȡһ��,��ŵ�p_line����
    if (fgets(p_line, MAXCLM, stdin) == NULL) return 0;//û�ж�ȡ�����ݣ��ļ��Ѿ�����������0

    //Ѱ�һ��з������û�ҵ����з�����ôҪô�ļ�������Ҫô��ʽ��������
    if ((p = strchr(p_line, '\n')) == NULL) {
        
        //�������ļ�������������0
        if (*p_line == '\032') return 0;           

        //�������û�л��з���ͬʱ���ǿ��У������ļ���δ��������˸��и�ʽ����ȷ��������һ��epoch line
        if (*p_line != '\0' || feof(stdin) == 0) {
            //������������У����ò��ܹ�������ʽ����ȷ���У���δ�Ӳ�����-s�������˳�����
            if (!skip_strange_epoch) error_exit(12, p_line);

            //������һ��epoch line
            skip_to_next(p_line);

            //���ض����˴����һ��
            return 2;
        }
        //������в��ǿ���ͬʱ�ļ����Ѿ��������ļ�ȱ�ٽ�����ʶ����ʾȱ�ٽ�����ʶ����ȷ�����ļ�����
        else {
            fprintf(stderr, "WARNING: null characters are detected at the end of file --> neglected.\n");
            exit_status = EXIT_WARNING;
            return 0;
        }
    }


    //���ڻ��з�

    //������\r\n�滻Ϊ\n
    if (*(p - 1) == '\r') { *(--p) = '\0'; };

    //�ӻ��з�\n�����е����һ���ַ�����ʼ��ǰ�ң����ҵ��ĵ�һ�����ǿո��ַ���ߵĵ�һ���ո��ó�\0
    while (*--p == ' ' && p > p_line) {}; *++p = '\0';        

    if (rinex_version == 2) {
        if (strlen(p_line) < 29 || *p_line != ' '
            || *(p_line + 27) != ' ' || !isdigit(*(p_line + 28))
            || (*(p_line + 29) != ' ' && !isdigit(*(p_line + 29)) && *(p_line + 29) != '\0')) {
            //���и�ʽ����ȷ��ͬʱδ��������-s��,�˳�
            if (!skip_strange_epoch) error_exit(6, p_line);

            //���и�ʽ����ȷ���жϵ�19���ַ��Ƿ�Ϊ��.��,���������������ڴ�
            if (*(p_line + 18) != '.')  CLEAR_BUFF;
            //��ת����һ��epoch line
            skip_to_next(p_line);
            return 2;
        }
    }
    //�汾Ϊ3��4ʱ
    else {   
        //�����һ���ַ����ǡ�>����������һ��epoch line
        if (*p_line != '>') {
            if (!skip_strange_epoch) error_exit(6, p_line);
            CLEAR_BUFF;
            skip_to_next(p_line);
            return 2;
        }

        //�����Ȳ��䵽42����43λ�óɡ�\0��
        while (p < (p_line + 41)) *p++ = ' '; /*** pad blank ***/
        *p = '\0';
    }
    return 1;
}
//������һ��epoch line
void skip_to_next(char* p_line) {
    //�����һ������ȷ��ʽ��
    fprintf(stderr, " WARNING at line %ld: strange format. skip to next epoch.\n", nl_count);
    exit_status = EXIT_WARNING;

    //Ѱ����һ��epoch line
    if (rinex_version == 2) {
        //����汾��Ϊ2��Ҫ���ݸ��ж�����ݸ�ʽ�Ƿ��ڷ�Χ�ڻ���ֵ�Ƿ���ϸ�ʽ�����жϸ����Ƿ�Ϊ�µ�epoch line
        do {                              
            read_chk_line(p_line);
        } while (strlen(p_line) < 29 || *p_line != ' ' || *(p_line + 3) != ' '
            || *(p_line + 6) != ' ' || *(p_line + 9) != ' '
            || *(p_line + 12) != ' ' || *(p_line + 15) != ' '
            || *(p_line + 26) != ' ' || *(p_line + 27) != ' '
            || !isdigit(*(p_line + 28)) || !isspace(*(p_line + 29))
            || (strlen(p_line) > 68 && *(p_line + 70) != '.'));
    }
    else {    
        //����汾�Ų���2���µ�epoch�п�ʼ�ַ�Ϊ��>��
        do {
            read_chk_line(p_line);
        } while (*p_line != '>');
    }
    //���������³�ʼ��
    initialize_all(oldline, &nsat_old, 0);             
}
//��ʼ��һЩ����
void initialize_all(char* oldline, int* nsat_old, int count) {
    strcpy(oldline, "&");        /**** initialize the epoch data arc ****/
    clk_order = -1;             /**** initialize the clock data arc ****/
    *nsat_old = 0;              /**** initialize the all satellite arcs ****/
    ep_count = count;
}
/*---------------------------------------------------------------------*/
void put_event_data(char* p_line) {
    /**** This routine is called when event flag >1 is set.  ****/
    /****      read # of event information lines and output  ****/
    int i, n;
    char* p;

    if (rinex_version == 2) {
        if (*(p_line + 26) == '.') error_exit(6, p_line);
        printf("&%s\n", (p_line + 1));
        if (strlen(p_line) > 29) {
            n = atoi((p_line + 29));     /** n: number of lines to follow **/
            for (i = 0; i < n; i++) {
                read_chk_line(p_line);
                printf("%s\n", p_line);
                if (strncmp((p_line + 60), "# / TYPES OF OBSERV", C1 * 19) == 0 && *(p_line + 5) != ' ') {
                    *flag[0] = '\0';
                    ntype = atoi(p_line);
                    if (ntype > MAXTYPE) error_exit(16, p_line);
                }
            }
        }
    }
    else {
        if (strlen(p_line) < 35 || *(p_line + 29) == '.') error_exit(6, p_line);
        /* chop blanks that were padded in get_next_epoch */
        p = strchr(p_line + 35, '\0'); while (*--p == ' ') {}; *++p = '\0';
        printf("%s\n", p_line);
        n = atoi((p_line + 32));         /** n: number of lines to follow **/
        for (i = 0; i < n; i++) {
            read_chk_line(p_line);
            printf("%s\n", p_line);
            if (strncmp((p_line + 60), "SYS / # / OBS TYPES", C1 * 19) == 0 && *p_line != ' ') {
                *flag[0] = '\0';
                ntype_gnss[(unsigned int)*p_line] = atoi((p_line + 3));
                if (ntype_gnss[(unsigned int)*p_line] > MAXTYPE) error_exit(16, p_line);
            }
        }
    }
}
//��ȡʱ������,û��ʾ������
void read_clock(char* p_clock, int shift_clk) {
    /****  read the clock offset value ****/
    /**  *p_clock : pointer to beginning of clock data **/

    char* p_dot; 
    p_dot = p_clock + 2;
    if (*p_dot != '.')error_exit(7, p_clock);

    strncpy(p_dot, p_dot + 1, C1 * shift_clk);  /**** shift digits because of too  ****/
    *(p_dot + shift_clk) = '.';             /**** many digits for fractional part ****/
    sscanf(p_clock, "%ld.%ld", &clk1.u[0], &clk1.l[0]);
    if (*p_clock == '-' || *(p_clock + 1) == '-') clk1.l[0] = -clk1.l[0];
    if (clk_order < ARC_ORDER) clk_order++;
    *p_clock = '\0';
}
/*---------------------------------------------------------------------*/
void process_clock(void) {
    int i;
    for (i = 0; i < clk_order; i++) {
        clk1.u[i + 1] = clk1.u[i] - clk0.u[i];
        clk1.l[i + 1] = clk1.l[i] - clk0.l[i];
    }
}
//�����Աȱ��к�ǰһ�е�����˳��ǰһ�С�C01C02C03�������С�C02C01C03��
int  set_sat_table(char* p_new, char* p_old, int nsat_old, int* sattbl) {                                
    int i, j;
    char* ps;

    for (i = 0; i < nsat; i++, p_new += 3) {
        *sattbl = -1;
        ps = p_old;
        for (j = 0; j < nsat_old; j++, ps += 3) {
            if (strncmp(p_new, ps, C3) == 0) {
                *sattbl = j;
                break;
            }
        }
        /*** check double entry ***/
        for (j = i + 1, ps = p_new + 3; j < nsat; j++, ps += 3) {
            if (strncmp(p_new, ps, C3) == 0) {
                if (!skip_strange_epoch) error_exit(13, p_new);
                fprintf(stderr, "WARNING:Duplicated satellite in one epoch at line %ld. ... skip\n", nl_count);
                return 1;
            }
        }
        sattbl++;
    }
    return 0;
}
/*---------------------------------------------------------------------*/
int  read_more_sat(int n, char* p) {
    /**** read continuation line of satellite list (for RINEX2) ****/
    char line[MAXCLM];

    do {
        p += 36;
        if (read_chk_line(line)) return 1;
        /**** append satellite table ****/
        if (line[2] == ' ') {
            sprintf(p, "%s", &line[32]);
        }
        else {                        /*** for the files before clarification of format ***/
            sprintf(p, "%s", &line[0]); /*** by W.Gurtner (IGS mail #1577)                ***/
        }
        n -= 12;
    } while (n > 12);
    return 0;
}
/*---------------------------------------------------------------------*/
void data(int* sattbl) {
    /********************************************************************/
    /*  Function : output the 3rd order difference of data              */
    /*       u : upper X digits of the data                             */
    /*       l : lower 5 digits of the data                             */
    /*            ( y = u*100 + l/1000 )                                */
    /*   py->u : upper digits of the 3rd order difference of the data   */
    /*   py->l : lower digits of the 3rd order difference of the data   */
    /********************************************************************/
    data_format* py1;
    int  i, j, * i0;
    char* p;

    for (i = 0, i0 = sattbl; i < nsat; i++, i0++) {
        for (j = 0, py1 = dy1[i]; j < ntype_record[i]; j++, py1++) {
            if (py1->order >= 0) {       /*** if the numerical data field is non-blank ***/
                if (*i0 < 0 || dy0[*i0][j].order == -1) {
                    /**** initialize the data arc ****/
                    py1->order = 0; p_buff += sprintf(p_buff, "%d&", ARC_ORDER);
                }
                else {
                    take_diff(py1, &(dy0[*i0][j]));
                    if (labs(py1->u[py1->order]) > 100000) {
                        /**** initialization of the arc for large cycle slip  ****/
                        py1->order = 0; p_buff += sprintf(p_buff, "%d&", ARC_ORDER);
                    }
                }
                putdiff(py1->u[py1->order], py1->l[py1->order]);
            }
            else if (*i0 >= 0 && rinex_version == 2) {
                /**** CRINEX1 (RINEX2) initialize flags for blank field, not put '&' ****/
                flag0[*i0][j * 2] = flag0[*i0][j * 2 + 1] = ' ';
            }
            if (j < ntype_record[i] - 1) *p_buff++ = ' ';   /** ' ' :field separator **/
        }
        *(p_buff++) = ' ';  /* write field separator */
        if (*i0 < 0) {             /* if new satellite initialize all LLI & SN flags */
            if (rinex_version == 2) {
                p_buff = strdiff("", flag[i], p_buff);
            }
            else {          /*  replace space with '&' for CRINEX3(RINEX3)  */
                for (p = flag[i]; *p != '\0'; p++) *p_buff++ = (*p == ' ') ? '&' : *p;
                *p_buff++ = '\n'; *p_buff = '\0';
            }
        }
        else {
            p_buff = strdiff(flag0[*i0], flag[i], p_buff);
        }
    }
}
//����ַ��Ա�s1��s2�����s2��s1��Ӧ�ַ���ͬ��ds��ͬλ���ÿո������ͬ��s2λ����Ϊ�ո�����&������
//���������ds�б���s2��Ӧ�ַ�
//����:s1Ϊ��2023 12 30 12 12 00�� s2Ϊ��2023 12 30  9 12 30��ת���� dsΪ��           &9    3 ��
char* strdiff(const char* s1, char* s2, char* ds) {
    //�Ա�s2��s1�����s2��s1��Ӧ�ַ���ͬ��ds��ͬλ���ÿո������ͬ��s2λ����Ϊ�ո�����&������
    for (; *s1 != '\0' && *s2 != '\0'; s2++) {
        if (*s2 == *(s1++))
            *ds++ = ' ';
        else if (*s2 == ' ')
            *ds++ = '&';
        else
            *ds++ = *s2;
    }

    //���s1��s2������s1ʣ�µķŵ�ds�У�Ȼ����Щ�ַ����κβ�Ϊ�ո���ַ���Ϊ&
    //����s1Ϊ��2023 12 30 12 12 00 54����s2Ϊ��2023 12 30  9 12 30����ת�����dsΪ��           &9    3  &&��
    strcpy(ds, s1);
    for (; *ds; ds++) { if (*ds != ' ') *ds = '&'; }

    //���s2��s1������s2ʣ�µ�ֱ��ȫ���ŵ�ds��
    while (*s2) *ds++ = *s2++;

    //���ds��󼸸��ַ��ǿո�����ʡ��
    for (ds--; *ds == ' '; ds--);    /*** find pointer of last non-space character ***/
    *++ds = '\n'; *++ds = '\0';       /*** chop spaces at the end of the line ***/
    return ds;
}
/*---------------------------------------------------------------------*/
int  ggetline(data_format* py1, char* flag, char* sat_id, int* ntype_rec) {
    /**** read data line for one satellite and       ****/
    /**** set data difference and flags to variables ****/ 
   
    char line[MAXCLM], * p, * pmax, * p_1st_rec;
    int i, j, nfield, max_field;

    if (read_chk_line(line)) return 1;

    //�����ļ�ͷ�е����ã�ȷ�������������ռ�ó���
    if (rinex_version == 2) {             /** for RINEX2 **/
        max_field = 5;                             /** maximum data types in one line **/
        *ntype_rec = ntype;                        /** # of data types for the satellite **/
        p_1st_rec = line;                          /** pointer to the start of the first record **/
    }
    else {                                /** for RINEX3 **/
        strncpy(sat_id, line, C3);                   /** put satellite ID to the list of satellites **/
        max_field = *ntype_rec = ntype_gnss[(unsigned int)line[0]];  /*** # of data types for the GNSS system ***/
        if (max_field < 0) {
            if (!skip_strange_epoch) error_exit(21, line);
            fprintf(stderr, "WARNING at line %ld. : GNSS type '%c' is not defined in the header. ... skip\n", nl_count, (unsigned int)line[0]);
            return 1;
        }
        p_1st_rec = line + 3;
    }
    for (i = 0; i < *ntype_rec; i += max_field) {                 
        
        //�ҵ��������ݵ����λ��
        nfield = (*ntype_rec - i < max_field ? *ntype_rec - i : max_field); 
        pmax = p_1st_rec + 16 * nfield;

        //�ҵ����н�����λ��
        p = strchr(line, '\0');
        //����������������λ��ǰ���Ѿ��������ÿո�������䵽pmax
        if (p < pmax) {
            while (p < pmax) *p++ = ' '; *p = '\0';
        }
        //������ݳ��ȴ���pmax�����������ݴ�������
        else {
            for (*p = ' '; p > pmax && *p == ' '; p--) {};
            if (p > pmax) {
                if (!skip_strange_epoch) error_exit(9, line);
                fprintf(stderr, "WARNING: mismatch of number of the data types at line %ld. ... skip\n", nl_count);
                return 1;
            }
        }

        
        for (j = 0, p = p_1st_rec; j < nfield; j++, p += 16, py1++) {
            //�����11��λ����Ϊ��.��֤���˴�������
            if (*(p + 10) == '.') {
                *flag++ = *(p + 14);
                *flag++ = *(p + 15);
                *(p + 14) = '\0';
                read_value(p, &(py1->u[0]), &(py1->l[0]));
                py1->order = 0;
            }
            //�ж��Ƿ�Ϊ�հ�����
            else if (strncmp(p, "              ", C14) == 0) {
                if (rinex_version == 2 && strncmp((p + 14), "  ", C2) != 0) error_exit(20, line);
                *flag++ = *(p + 14);
                *flag++ = *(p + 15);
                py1->order = -1;
            }
            //�����ļ���ʽ���⣬��Ҫ����
            else {
                if (!skip_strange_epoch) error_exit(10, p);
                fprintf(stderr, "WARNING: abnormal data field at line %ld....skip\n", nl_count);
                return 1;
            }
        }
        //��ȡ��һ��
        if (i + max_field < *ntype_rec) {
            if (read_chk_line(line)) return 1;   /* read continuation line */
        }
    }
    *flag = '\0';
    return 0;
}
//��ȡ���ݣ�������Ϊ���������ִ洢
void read_value(char* p, long* pu, long* pl) {
    /**** divide the data into lower 5 digits and upper digits     ****/
    /**** input p :  pointer to one record (14 characters + '\0')  ****/
    /**** output  *pu, *pl: upper and lower digits the data        ****/

    char* p7, * p8, * p9;
    p7 = p + 7;
    p8 = p7 + 1;
    p9 = p8 + 1;

    *(p9 + 1) = *p9;            /* shift two digits: ex. 123.456 -> 1223456,  -.345 ->   -345 */
    *p9 = *p8;                /*                       -12.345 -> -112345, -1.234 -> --1234 */
    *pl = atol(p9);           /*                         0.123 ->  . 0123, -0.123 -> --0123 */

    if (*p7 == ' ') {
        *pu = 0;
    }
    else if (*p7 == '-') {
        *pu = 0;
        *pl = -*pl;
    }
    else {
        *p8 = '.';
        *pu = atol(p);
        if (*pu < 0) *pl = -*pl;
    }
}
/*---------------------------------------------------------------------*/
void take_diff(data_format* py1, data_format* py0) {
    int k;

    py1->order = py0->order;
    if (py1->order < ARC_ORDER) (py1->order)++;
    if (py1->order > 0) {
        for (k = 0; k < py1->order; k++) {
            py1->u[k + 1] = py1->u[k] - py0->u[k];
            py1->l[k + 1] = py1->l[k] - py0->l[k];
        }
    }
}
/*---------------------------------------------------------------------*/
void putdiff(long dddu, long dddl) {

    dddu += dddl / 100000; dddl %= 100000;
    if (dddu < 0 && dddl>0) {
        dddu++; dddl -= 100000;
    }
    else if (dddu > 0 && dddl < 0) {
        dddu--; dddl += 100000;
    }

    if (dddu == 0) {
        p_buff += sprintf(p_buff, "%ld", dddl);
    }
    else {
        p_buff += sprintf(p_buff, "%ld%5.5ld", dddu, labs(dddl));
    }
}
/*---------------------------------------------------------------------*/
void put_clock(long du, long dl, int c_order) {
    /***********************************/
    /****  output clock diff. data  ****/
    /***********************************/
    du += dl / 100000000; dl %= 100000000;
    if (du < 0 && dl>0) {
        du++; dl -= 100000000;
    }
    else if (du > 0 && dl < 0) {
        du--; dl += 100000000;
    }
    if (c_order == 0) p_buff += sprintf(p_buff, "%d&", ARC_ORDER);
    if (du == 0) {
        p_buff += sprintf(p_buff, "%ld\n", dl);
    }
    else {
        p_buff += sprintf(p_buff, "%ld%8.8ld\n", du, labs(dl));
    }
}
//��һ�����ݣ������еĸ�ʽ�����е����һ���ַ�����Ϊ���з�������\n��
int  read_chk_line(char* line) {
    char* p;
    nl_count++;
    if (fgets(line, MAXCLM, stdin) == NULL) error_exit(11, line);
    if ((p = strchr(line, '\n')) == NULL) {
        if (fgetc(stdin) == EOF) {
            error_exit(11, line);
        }
        else {
            if (!skip_strange_epoch) error_exit(12, line);
            fprintf(stderr, "WARNING: null character is found or the line is too long (>%d) at line %ld.\n", MAXCLM, nl_count);
            return 1;
        }
    }
    if (*(p - 1) == '\r')p--;   /*** check DOS CR/LF ***/

    //������з�ǰ�����пո�����"aaaaaa aaa   a   \n"�����aaaaaa aaa   a\n��
    while (*--p == ' ' && p > line) {}; *++p = '\0'; 
    return 0;
}
/*---------------------------------------------------------------------*/
void error_exit(int error_no, char* string) {
    if (error_no == 1) {
        fprintf(stderr, "Usage: %s [file] [-] [-f] [-e # of epochs] [-s] [-d] [-h]\n", string);
        fprintf(stderr, "    stdin and stdout are used if input file name is not given.\n");
        fprintf(stderr, "    -       : output to stdout\n");
        fprintf(stderr, "    -f      : force overwrite of output file\n");
        fprintf(stderr, "    -e #    : initialize the compression operation at every # epochs\n");
        fprintf(stderr, "              When some part of the Compact RINEX file is lost, the data\n");
        fprintf(stderr, "              can not be recovered thereafter until all the data arc are\n");
        fprintf(stderr, "              initialized for differential operation. This option may be used to\n");
        fprintf(stderr, "              increase chances to recover parts of data by using an option of\n");
        fprintf(stderr, "              CRX2RNX(ver. 4.0 or after) with cost of increase of file size.\n");
        fprintf(stderr, "    -s      : warn and skip strange epochs (default: stop with error status)\n");
        fprintf(stderr, "    -d      : delete the input file if conversion finishes without errors\n");
        fprintf(stderr, "              (i.e. exit code = %d or %d).\n", EXIT_SUCCESS, EXIT_WARNING);
        fprintf(stderr, "              This option does nothing if stdin is used for the input.\n");
        fprintf(stderr, "    -h      : display this message\n\n");
        fprintf(stderr, "    exit code = %d (success)\n", EXIT_SUCCESS);
        fprintf(stderr, "              = %d (error)\n", EXIT_FAILURE);
        fprintf(stderr, "              = %d (warning)\n", EXIT_WARNING);
        fprintf(stderr, "    [version : %s]\n", VERSION);
        exit(EXIT_FAILURE);
    }
    if (error_no == 4) {
        fprintf(stderr, "ERROR : invalid file name  %s\n", string);
        fprintf(stderr, "The extension of the input file name should be [.??o] or [.rnx].\n");
        fprintf(stderr, "To convert the files whose name is not fit to the above conventions,\n");
        fprintf(stderr, "use of this program as a filter is also possible. \n");
        fprintf(stderr, "    for example)  cat file.in | %s - > file.out\n", PROGNAME);
        exit(EXIT_FAILURE);
    }
    if (error_no == 5) {
        fprintf(stderr, "ERROR : can't open %s\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 6) {
        fprintf(stderr, "ERROR when reading line %ld.\n", nl_count);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 7) {
        fprintf(stderr, "ERROR at line %ld: invalid format for clock offset.\n", nl_count);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 8) {
        fprintf(stderr, "ERROR at line %ld : number of satellites exceed the maximum(%d).\n", nl_count, MAXSAT);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 9) {
        fprintf(stderr, "ERROR at line %ld : mismatch of number of the data types.\n", nl_count);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 10) {
        fprintf(stderr, "ERROR at line %ld : abnormal data field.\n", nl_count);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 11) {
        fprintf(stderr, "ERROR : The RINEX file seems to be truncated in the middle.\n");
        fprintf(stderr, "        The conversion is interrupted after reading line %ld :\n", nl_count);
        fprintf(stderr, "        start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 12) {
        fprintf(stderr, "ERROR at line %ld. : null character is found or the line is too long (>%d).\n", nl_count, MAXCLM);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 13) {
        fprintf(stderr, "ERROR at line %ld. : Duplicated satellite in one epoch.\n", nl_count);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 14) {
        fprintf(stderr, "ERROR at line %ld. : Length of file name exceed MAXCLM(%d).\n", nl_count, MAXCLM);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 15) {
        fprintf(stderr, "The first line is :\n%s\n\n", string);
        fprintf(stderr, "ERROR : The file format is not valid. This program is applicable\n");
        fprintf(stderr, "        only to RINEX Version 2/3/4 Observation file.\n");
        exit(EXIT_FAILURE);
    }
    if (error_no == 16) {
        fprintf(stderr, "ERROR at line %ld. : Number of data types exceed MAXTYPE(%d).\n", nl_count, MAXTYPE);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 20) {
        fprintf(stderr, "ERROR at line %ld. : data is blank but there is flag.\n", nl_count);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
    if (error_no == 21) {
        fprintf(stderr, "ERROR at line %ld. : GNSS type '%c' is not defined in the header.\n", nl_count, (unsigned int)string[0]);
        fprintf(stderr, "     start>%s<end\n", string);
        exit(EXIT_FAILURE);
    }
}
/*---------------------------------------------------------------------*/
void no_error_exit() {
    if (delete_if_no_error && exit_status != 2 && n_infile == 1)  remove(infile);
    exit(exit_status);
}
