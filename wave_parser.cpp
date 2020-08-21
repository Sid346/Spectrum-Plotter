#include "mainwindow.h"


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
    int data_size;
}head;

QVector<double>* data_read(QString filename, int* freq) {
    QFile file(filename);

    file.open(QIODevice::ReadOnly);
    head *header = new head;
    char strm;
    file.read((char*)header,sizeof(head));

    *freq  =  header->sample_rate;
    int i = header->data_size;
    QVector<double>* data = new  QVector<double>;

    while (i >0)
    {
        file.read(&strm, sizeof(short));
        short a = (qFromLittleEndian<short>((uchar*)&strm));
        *data << (a/32767.00);
        i -= 2;
    }
    file.close();
    return data;
}
