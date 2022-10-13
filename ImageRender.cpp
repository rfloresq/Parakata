// ******************************************************************
// Purpose: Definition of ImageRender class functions.
//
// Roberto Flores-Moreno
// ******************************************************************
//
#include <QtGui>

#include <QtOpenGL/QGLWidget>
#include <ImageRender.h>
#include <PlotManager.h>
#include <Parakata.h>
#include <Viewer.h>
#include <global.h>

#define AGIF_ROT 0

static void system_image( QPixmap *pixmap , char* format, char* filename )
{
  char command[1024];

  if ( QString( format ) == "png" ||
       QString( format ) == "bmp" )
  {
    pixmap->save( filename , format , 0 );
  }
  else
  {
// KPU 2011
    pixmap->save( "imagen.png" , "png" , 100 );
  //  sprintf( command , "convert /tmp/tmp.png /tmp/tmp.%s", format );
    system( command );
    system("rm -f /tmp/tmp.png" );
    sprintf( command , "mv /tmp/tmp.%s %s", format , filename );
    system( command );
  };
}

ImageRender::ImageRender( Viewer *iviewer )
           : QWidget( 0 ) 
{
        viewer = iviewer;
        QVBoxLayout *mainLayout = new QVBoxLayout;
//
        QGroupBox *labelGroupBox = new QGroupBox(tr("Label"));
        QHBoxLayout *labellayout = new QHBoxLayout;
        label = new QLineEdit( labelGroupBox ); 
        label->setText( QString( "" ) ); 
       // connect( label , SIGNAL( returnPressed() ),
       //          this , SLOT( RedrawGLWidget() ));
        labellayout->addWidget( label );
        labelGroupBox->setLayout( labellayout );
        mainLayout->addWidget( labelGroupBox );
//
        QGroupBox *horizontalGroupBox = new QGroupBox(tr("Axis"));
        QHBoxLayout *axislayout = new QHBoxLayout;
        for (int i = 0; i < 3; i++) 
        {
          axised[i] = new QLineEdit( horizontalGroupBox ); 
          axised[i]->setText( QString::number( 1.0 ) ); 
          axislayout->addWidget( axised[i] );
        };
        horizontalGroupBox->setLayout( axislayout );
        mainLayout->addWidget( horizontalGroupBox );
//
        QHBoxLayout *descLayout = new QHBoxLayout;
        formatCombo = new QComboBox;
        formatCombo->addItem( QString( "Animated GIF" ));
        formatCombo->addItem( QString( "PS" ));
        formatCombo->addItem( QString( "GIF" ));
        formatCombo->addItem( QString( "PDF" ));
        formatCombo->addItem( QString( "PNG" ));
        formatCombo->addItem( QString( "JPG" ));
        formatCombo->addItem( QString( "JPEG" ));
        formatCombo->addItem( QString( "BMP" ));
        descLayout->addWidget( formatCombo );
        connect( formatCombo, SIGNAL(activated(const QString &)),
                 this, SLOT( ChangeFormat(const QString &)));
//
        QPushButton *writeB = new QPushButton( tr("Write") );
        connect( writeB, SIGNAL(clicked()), this, SLOT(Write()));
        descLayout->addWidget( writeB );
//
        QPushButton *doneB = new QPushButton( tr("Done") );
        connect( doneB, SIGNAL(clicked()), this, SLOT(close()));
        descLayout->addWidget( doneB );
//
        mainLayout->addLayout( descLayout );
        setLayout( mainLayout );
//
        setWindowTitle( "Image Renderer" );

}

ImageRender::~ImageRender()
{
}

bool ImageRender::GetFileName( char* format )
{
    //    QString initialPath = QDir::currentPath() + "/untitled." + format;
       QString initialPath = QDir::currentPath();
        QString fn = QFileDialog::getSaveFileName(this, 
        "Save Image to File", initialPath,format);
        if ( fn.isEmpty() || fn.isNull() ) 
        {
          return false;
        }
        else
        {
          sprintf(filename,"%s",(fn.toStdString()).c_str());
          return true;
        }
}

void ImageRender::DoAnimatedGIF( int type )
{
        if ( GetFileName( "gif" ) )
        {
//
          if ( type == AGIF_ROT )
          {
            int n = 10;
            char imgname[1024];
            char command[1024];
//
            prk->Report("Writing GIF image file %s ...", filename );
            for ( int i = 0 ; i < n ; i++ )
            {
              sprintf(imgname,"/tmp/%simg%d.gif",PROGRAM_NAME,i);
              viewer->Rotate( 36.0 , 
                              (axised[0]->text()).toDouble() , 
                              (axised[1]->text()).toDouble() , 
                              (axised[2]->text()).toDouble() );
              viewer->Redraw();
              QPixmap pixmap = viewer->renderPixmap( viewer->width(), 
                                              viewer->height());
              system_image( &pixmap , "gif" , imgname );
              prk->Report("Writing GIF image file %s ...(%d of %d)",
                                    filename , i + 1 , n );
            };
            sprintf(command,"convert -delay 5 /tmp/%simg*.gif %s",
                    PROGRAM_NAME,filename);
            system( command );
            sprintf(command,"mv /tmp/%simg*.gif /work/",PROGRAM_NAME);
            system( command );
            prk->Report("Writing GIF image file %s ... DONE",
                                  filename );
          };
        };
}

void ImageRender::ChangeFormat( const QString& newfmt )
{
        if ( newfmt == "Animated GIF" )
        {
          axised[0]->setEnabled( true );
          axised[1]->setEnabled( true );
          axised[2]->setEnabled( true );
        }
        else
        {
          axised[0]->setEnabled( false );
          axised[1]->setEnabled( false );
          axised[2]->setEnabled( false );
        };
}

void ImageRender::Write( )
{
        viewer->pmgr->SetUsageOfLists( false );
        QString text = formatCombo->currentText();
        if ( text == "Animated GIF" )
        {
          DoAnimatedGIF( AGIF_ROT );
        }
        else
        {
          char format[64];
          sprintf( format , "%s", ((text.toLower()).toStdString()).c_str() );
          if ( GetFileName( format ) )
          {
            prk->Report( "Writing %s image file %s ...",
                                   format ,  filename );
            QPixmap pixmap = viewer->renderPixmap(viewer->width(), 
                                                  viewer->height());
            system_image( &pixmap , format , filename );
            prk->Report( "Writing %s image file %s ... DONE",
                                   format , filename);
          };
        };
        viewer->pmgr->SetUsageOfLists( true );
}

