#include "sd_card.h"
#include "main.h"
#include "fatfs.h"

char TxBuffer[250];
FRESULT FR_Status;


FATFS *FS_Ptr;
DWORD FreeClusters;
uint32_t TotalSize, FreeSpace;

void mount_sd(){
	FATFS FatFs;

	FR_Status = f_mount(&FatFs, "", 1);
	if (FR_Status != FR_OK)
	{
	  sprintf(TxBuffer, "Error! While Mounting SD Card, Error Code: (%i)\r\n", FR_Status);
	  return;
	}
	sprintf(TxBuffer, "SD Card Mounted Successfully! \r\n\n");
}

void write_sd(char* file_name, char* data){
	FIL Fil;

	//Open the file
	FR_Status = f_open(&Fil, file_name, FA_WRITE | FA_OPEN_ALWAYS);
	if(FR_Status != FR_OK)
	{
	  sprintf(TxBuffer, "Error! While Creating/Opening Text File, Error Code: (%i)\r\n", FR_Status);
	  return;
	}
	FR_Status = f_lseek(&Fil, f_size(&Fil)); // Move The File Pointer To The EOF (End-Of-File)
	if(FR_Status != FR_OK)
	{
	  sprintf(TxBuffer, "Error! While moving pointer to end of file, Error Code: (%i)\r\n", FR_Status);
	  return;
	}

	sprintf(TxBuffer, "Text File Opened! Writing Data To The Text File..\r\n\n");
	// (1) Write Data To The Text File [ Using f_puts() Function ]
	f_puts(data, &Fil);
	// Close The File
	f_close(&Fil);
}

void read_sd(char* file_name, uint32_t line_number, char* RW_buffer, size_t buffer_size)
{
    FIL Fil;
    UINT bytesRead;
    uint32_t current_line = 0;

    // Open the file for reading
    FR_Status = f_open(&Fil, file_name, FA_READ);
    if(FR_Status != FR_OK) {
        sprintf(TxBuffer, "Error opening file %s, code: %i\r\n", file_name, FR_Status);
        RW_buffer[0] = '\0';
        return;
    }

    RW_buffer[0] = '\0';
    size_t pos = 0;

    while(f_gets(&RW_buffer[pos], buffer_size - pos, &Fil)) {
        if(current_line == line_number) {
            // Remove newline characters
            size_t len = strlen(RW_buffer);
            if(len > 0 && (RW_buffer[len-1] == '\n' || RW_buffer[len-1] == '\r')) {
                RW_buffer[len-1] = '\0';
                if(len > 1 && RW_buffer[len-2] == '\r') RW_buffer[len-2] = '\0';
            }
            break;
        }

        // Move to next line
        current_line++;
    }

    // If requested line does not exist, return empty string
    if(current_line < line_number) {
        RW_buffer[0] = '\0';
    }

    f_close(&Fil);
}

void delete_file(char* file_name){
	FR_Status = f_unlink(file_name);
	if (FR_Status != FR_OK){
		sprintf(TxBuffer, "Error! While Deleting The (TextFileWrite.txt) File.. \r\n");
	}
}

void unmount_sd(char* file_name){
	FR_Status = f_mount(NULL, "", 0);
	if (FR_Status != FR_OK)
	{
		sprintf(TxBuffer, "Error! While Un-mounting SD Card, Error Code: (%i)\r\n", FR_Status);
	} else{
		sprintf(TxBuffer, "SD Card Un-mounted Successfully! \r\n");
	}
}
