#include "test.h"


struct header
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
    char data_id[4];    					// "data" written in ascii
}head;

QVector<double>* data_read(QString filename) {
    QFile file(filename);
    QVector<double> *data = new QVector<double>;
    file.open(QIODevice::ReadOnly);
    QByteArray line = file.read(sizeof head);
    qDebug() << line << endl;
    char strm;
    file.read(&strm,4);
    //qDebug() << qFromLittleEndian<quint32>((uchar*)&strm);

    while (!file.atEnd())
    {
        file.read(&strm,2);
        if (qFromLittleEndian<short>((uchar*)&strm))
            *data << (qFromLittleEndian<short>((uchar*)&strm))/65535.0;
    }
    file.close();
    return data;
}
