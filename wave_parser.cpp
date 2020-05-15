#include "test.h"


typedef struct header
{
    char fileformat[4];
    int file_size;
    char subformat[4];
    char subformat_id[4];
    int chunk_bits;     					// 16or18or40 due to pcm it is 16 here
    short int audio_format;    					// little or big endian
    short int num_channels;     				// 2 here for left and right
    int sample_rate;						// sample_rate denotes the sampling rate.
    int byte_rate;           					// bytes  per second
    short int bytes_per_frame;
    short int bits_per_sample;
    char data_id[4];
}head;

QVector<double>* data_read(QString filename, int* freq) {
    QFile file(filename);
    QVector<double> *data = new QVector<double>;
    
    file.open(QIODevice::ReadOnly);
    head *header;
    char strm;
    file.read((char*)header,sizeof head);
  
    *freq  =  header->sample_rate;
    file.read(&strm, sizeof(short));
    while (!file.atEnd())
    {
        
        file.read(&strm, sizeof(short));
        if (qFromLittleEndian<short>((uchar*)&strm)) {
            short a = (qFromLittleEndian<short>((uchar*)&strm));
            *data << (a/32778.00);
           
        }
    }
    file.close();
    return data;
}
