// ******************************************************************
// Purpose: Viewer object constructor/destructor.
//
// 2003-2021, Roberto Flores-Moreno
// ******************************************************************

#include <math.h>
#include <macros.h>

#include<iostream>

#include <QtOpenGL/QtOpenGL>
#include <QWheelEvent>

#include <fcns.h>
#include <global.h>

#include <Parakata.h>
#include <Viewer.h>
#include <GeometryEditor.h>
#include <ElementSelector.h>
#include <PlotManager.h>
#include <ImageRender.h>
#include <Vector.h>

bool perspective = false;

Viewer::Viewer( Parakata *is )
      : QGLWidget( is )
{
  prk = is;
  image_render = new ImageRender( this );
  pmgr = new PlotManager( prk , this );

  ResetRotation();

  trans[0] = 0.0; 
  trans[1] = 0.0; 
  trans[2] = -5.0; 
  scale = 1.0;

  bg[0] = 1.0;
  bg[1] = 1.0;
  bg[2] = 1.0;
  bg[4] = 0.0;

  tag_type = "None";
  tag_color[0] = 0.7;
  tag_color[1] = 0.7;
  tag_color[2] = 0.7;
  tag_color[3] = 0.0;

  arrow_type = "None";
  arrow_color[0] = 1.0;
  arrow_color[1] = 0.0;
  arrow_color[2] = 0.0;
  arrow_color[3] = 0.0;

  cp_color[0] = 0.0;
  cp_color[1] = 1.0;
  cp_color[2] = 0.0;

  view_width = width()/10;
  view_height = view_width*((float)height())/(float)width();
  view_near = -100.0;
  view_far = 100.0;
  perspective = false;
  if (perspective) 
  {
    view_near = -10.0;
    view_far = 30.0;
  }
  q = gluNewQuadric();

  this->setFocusPolicy(Qt::StrongFocus);

  user_label.clear();
  selected_atoms.clear();
  monitored_bonds.clear();
  monitored_angles.clear();
  monitored_dihedrals.clear();
}

Viewer::~Viewer()
{
  user_label.clear();
  selected_atoms.clear();
  monitored_bonds.clear();
  monitored_angles.clear();
  monitored_dihedrals.clear();
  delete pmgr;
}

void Viewer::Project()
{
  makeCurrent();
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  if (perspective) 
    glFrustum( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  else
    glOrtho( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
}
 
void Viewer::initializeGL( void )
{
  GLfloat whiteDir[4] = {2.0, 2.0, 2.0, 1.0};
  GLfloat lightPos0[4] = {30.0,30.0,30.0,1.0};

  glMaterialf(GL_FRONT, GL_SHININESS, 80.0);
  glLightfv( GL_LIGHT0 , GL_POSITION , lightPos0 );
  glEnable( GL_DEPTH_TEST );
  glEnable( GL_LIGHTING );
  glEnable( GL_LIGHT0 );
  glEnable( GL_COLOR );
  glEnable( GL_COLOR_MATERIAL );
  glShadeModel( GL_SMOOTH );
  glMaterialfv( GL_FRONT, GL_SPECULAR, whiteDir);
  glEnable(GL_NORMALIZE);
  glClearColor( bg[0] , bg[1] , bg[2] , bg[3] );

  resizeGL( width() , height() );
}

void Viewer::resizeGL( int w , int h )
{
  glViewport( 0, 0, w, h );
  view_height = view_width*((float)h/(float)w);
  Project();
}

void Viewer::Redraw()
{
  paintGL();
}

void Viewer::paintGL()
{
  Project();
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glMatrixMode( GL_MODELVIEW );
  glLoadIdentity();

  glTranslatef( trans[0] , trans[1] , trans[2] );
  glPushMatrix();
  glMultMatrixd( rot );
  glScalef( scale , scale , scale );
 
  DrawMolecule( 0 , NATOM-1 , geosrs , geocrs );
  //int half = NATOM/2;
  //DrawMolecule( 0 , half-1 , geosrs , geocrs );
  //for (int i=0;i<MAXEL;i++)
  //  for (int j=0;j<3;j++)
  //    ELEMENT_COLOR[i][j] = ELEMENT_COLOR[i][j]*0.3;
  //DrawMolecule( half , NATOM-1 , geosrs , geocrs );
  //for (int i=0;i<MAXEL;i++)
  //  for (int j=0;j<3;j++)
  //    ELEMENT_COLOR[i][j] = ELEMENT_COLOR[i][j]/0.3;

  DrawSurfaces();

  DrawStrings();

  glPopMatrix();

  swapBuffers();
}

void Viewer::mousePressEvent( QMouseEvent *e )
{
  mouse[0] = e->x();
  mouse[1] = e->y();
  if ( ( prk->status & SFB  ) ||
       ( prk->status & SFA  ) ||
       ( prk->status & SFD  ) ||
       ( prk->status & SFFV ) ||
       ( prk->status & SFRV ) ||
       ( prk->status & SFAV ) ||
       ( prk->status & SFDV ) ||
       ( prk->status & SFPB ) ||       
       ( prk->status & SFPG ) )
  {
    PickAtom();
  };
  if ( prk->status & SFSURF )
  {
    PickSurface();
  }
}

void Viewer::mouseMoveEvent( QMouseEvent *e )
{
  bool rotate = false;
  int sx = e->x() - mouse[0];
  int sy = e->y() - mouse[1];

  if ( abs(sx) > (signed)width()/100 && abs(sx) > abs(sy) ) 
  {
    mouse[0] = e->x();
    Rotate( (double)((180*sx)/(signed)width() % 360), 0.0, 1.0, 0.0); 
    rotate = true;
  };

  if ( abs(sy) > (signed)height()/100 && abs(sx) < abs(sy) ) 
  {
    mouse[1] = e->y();
    Rotate( (double)((180*sy)/(signed)height() % 360), 1.0, 0.0, 0.0); 
    rotate = true;
  };

  if ( rotate )
  {
    Redraw();
  }
}

void Viewer::mouseReleaseEvent( QMouseEvent *e )
{
  mouse[0] = e->x();
  mouse[1] = e->y();
}

void Viewer::keyPressEvent(QKeyEvent *keyEvent)
{
  if (keyEvent->key()==Qt::Key_Up)
  {
    Rotate( -10.0, 1.0, 0.0, 0.0); 
    Redraw();
    keyEvent->accept();
  }
  else if (keyEvent->key()==Qt::Key_Down)
  {
    Rotate( 10.0, 1.0, 0.0, 0.0); 
    Redraw();
    keyEvent->accept();
  }
  else if (keyEvent->key()==Qt::Key_Left)
  {
    Rotate( -10.0, 0.0, 1.0, 0.0); 
    Redraw();
    keyEvent->accept();
  }
  else if (keyEvent->key()==Qt::Key_Right)
  {
    Rotate( 10.0, 0.0, 1.0, 0.0); 
    Redraw();
    keyEvent->accept();
  }
  else
  {
    keyEvent->ignore();
  }
}

void Viewer::Rotate( double angle , double x , double y , double z )
{
  int i,j,k;
  double cosp,sinp,cost,norma;
  double m[3][3];
  double mm[3][3];
  double rm[3][3];
  double u[3];

  u[0] = x;
  u[1] = y;
  u[2] = z;
  norma = sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  if ( norma == 0.0 )
  {
    return;
  }
  u[0] = u[0]/norma;
  u[1] = u[1]/norma;
  u[2] = u[2]/norma;
  angle = angle*M_PI/180.0;
  cosp = cos(angle);
  sinp = sin(angle);
  cost = 1.0 - cosp;
  m[0][0] = u[0]*u[0]*cost + cosp;
  m[0][1] = u[0]*u[1]*cost - u[2]*sinp;
  m[0][2] = u[0]*u[2]*cost + u[1]*sinp;
  m[1][0] = u[0]*u[1]*cost + u[2]*sinp;
  m[1][1] = u[1]*u[1]*cost + cosp;
  m[1][2] = u[1]*u[2]*cost - u[0]*sinp;
  m[2][0] = u[0]*u[2]*cost - u[1]*sinp;
  m[2][1] = u[1]*u[2]*cost + u[0]*sinp;
  m[2][2] = u[2]*u[2]*cost + cosp;

  k = 0;
  for ( i = 0 ; i < 3 ; i++ )
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      rm[j][i] = rot[k];
      k++;
    };
    k++;
  };

  for ( i = 0 ; i < 3 ; i++ )
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      mm[i][j] = 0.0;
      for ( k = 0 ; k < 3 ; k++ )
      {
        mm[i][j] += m[i][k]*rm[k][j];
      };
    };
  };

  k = 0;
  for ( i = 0 ; i < 3 ; i++ )
  {
    for ( j = 0 ; j < 3 ; j++ )
    {
      rot[k] = mm[j][i];
      k++;
    };
    k++;
  };

}
void Viewer::ResetRotation(void)
{
  for ( int i = 0 ; i < 16 ; i++ )
    rot[i] = 0.0;
  rot[0]  = 1.0;
  rot[5]  = 1.0;
  rot[10] = 1.0;
  rot[15] = 1.0;
}

void Viewer::ZoomIn()
{
  scale *= 1.1;
  Redraw();
}

void Viewer::ZoomOut()
{
  scale /= 1.1;
  Redraw();
}

void Viewer::wheelEvent(QWheelEvent *event)
{
    if (event->delta()>0){
        ZoomIn();
    }else{
        if (event->delta()<0){
            ZoomOut();
        }
    }
}

void Viewer::PickAtom( void )
{
  GLfloat picksize=3.0;
  GLuint buffer[512];
  GLint viewport[4];
  GLint hits;
  int iatom;

  glGetIntegerv( GL_VIEWPORT , viewport );
  glSelectBuffer( 512 , buffer );
  glInitNames();
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPickMatrix( (GLfloat)mouse[0] ,
                 (GLfloat)(viewport[3]-mouse[1]) ,
                 picksize , picksize , viewport );
  if (perspective) 
    glFrustum( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  else
    glOrtho( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  glRenderMode( GL_SELECT );
  glMatrixMode( GL_MODELVIEW );
  glTranslatef( trans[0] , trans[1] , trans[2] );
  glPushMatrix();
  glMultMatrixd( rot );
  glScalef( scale , scale , scale );
  for ( iatom = 0 ; iatom < NATOM ; iatom++ )
  {
    glPushMatrix();
    glPushName( iatom );
    DrawMolecule( iatom , iatom , geosrs , 0.0 );
    glPopName();
    glPopMatrix();
  }

  glPopMatrix();
  hits = glRenderMode( GL_RENDER );
  Project();
  if ( hits > 0 )
  {
    SaveSelection( (int) buffer[3] );
  }; 
}

void Viewer::PickSurface( void )
{
  GLfloat picksize=3.0;
  GLuint buffer[512];
  GLint viewport[4];
  GLint hits;

  glGetIntegerv( GL_VIEWPORT , viewport );
  glSelectBuffer( 512 , buffer );
  glInitNames();
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  gluPickMatrix( (GLfloat)mouse[0] ,
                 (GLfloat)(viewport[3]-mouse[1]) ,
                  picksize , picksize , viewport );
  if (perspective) 
    glFrustum( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  else
    glOrtho( -view_width/2 , view_width/2, 
           -view_height/2, view_height/2,
            view_near, view_far );
  glRenderMode( GL_SELECT );
  glMatrixMode( GL_MODELVIEW );
  glTranslatef( trans[0] , trans[1] , trans[2] );
  glPushMatrix();
  glMultMatrixd( rot );
  glScalef( scale , scale , scale );

  for ( int i = 0 ; i < pmgr->NSurf() ; i++ )
  {
    glPushMatrix();
    glPushName( i );
    pmgr->PlotSurface( i );
    glPopName();
    glPopMatrix();
  };

  glPopMatrix();
  hits = glRenderMode( GL_RENDER );
  Project();
  if ( hits > 0 )
  {
    SaveSelection( (int) buffer[3] );
  }; 
}

void Viewer::DrawSphere( double* p , double radi , int n )
{
  glTranslatef ( p[0] , p[1] , p[2] );
  gluSphere( q , radi , n, n/2 );
  glTranslatef (-p[0] ,-p[1] ,-p[2] );
}

void Viewer::DrawCylinder(double *start, double *end, 
                          double radi , int n )
{
  GLdouble v[3];
  GLdouble dval, ca, sa, cb, sb;
  GLdouble tmat[4][4];
  GLfloat height; 

  v[0] = (GLdouble)(end[0]-start[0]); 
  v[1] = (GLdouble)(end[1]-start[1]); 
  v[2] = (GLdouble)(end[2]-start[2]); 
  height = (GLfloat) sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  dval = sqrt(v[0]*v[0] + v[1]*v[1]); 
  if (dval != 0.0) 
  {
    ca = v[1]/dval;
    sa = v[0]/dval;
    v[1] = v[0]*sa + v[1]*ca;
  } 
  else
  {
    ca = 1.0;
    sa = 0.0;
  };

  dval = sqrt(v[1]*v[1] + v[2]*v[2]); 
  if (dval != 0.0) 
  {
    cb = v[2]/dval;
    sb = v[1]/dval;
  }
  else 
  {
    cb = 1.0;
    sb = 0.0;
  };

  tmat[0][0] = ca;   
  tmat[1][0] = sa*cb; 
  tmat[2][0] = sa*sb; 
  tmat[3][0] = (GLdouble)start[0];
  tmat[0][1] = (-sa); 
  tmat[1][1] = ca*cb; 
  tmat[2][1] = ca*sb; 
  tmat[3][1] = (GLdouble)start[1];
  tmat[0][2] = 0.0;   
  tmat[1][2] = (-sb); 
  tmat[2][2] = cb;    
  tmat[3][2] = (GLdouble)start[2];
  tmat[0][3] = 0.0;   
  tmat[1][3] = 0.0;   
  tmat[2][3] = 0.0;   
  tmat[3][3] = 1.0;

  glPushMatrix();
  glMultMatrixd((GLdouble*)tmat);
  gluCylinder( q, (GLfloat)radi, (GLfloat)radi, height, n, 1);
  glPopMatrix();
}

void Viewer::DrawAtom( int atom , double radi , int n )
{
  double p[3];

  p[0] = COORD[1][atom][0];
  p[1] = COORD[1][atom][1];
  p[2] = COORD[1][atom][2];

  glColor3dv( ELEMENT_COLOR[ATOMNO[atom]] );

  DrawSphere( p , radi , n );

  if ( SelectedAtom( atom ) >= 0 )
  {
    double col[4] = { 0.0 , 0.0 ,0.0 , 0.0};
    col[0] = 1.0 - ELEMENT_COLOR[ATOMNO[atom]][0];
    col[1] = 1.0 - ELEMENT_COLOR[ATOMNO[atom]][1];
    col[2] = 1.0 - ELEMENT_COLOR[ATOMNO[atom]][2];
    glColor4dv( col );
    glLineWidth( 2 );
    gluQuadricDrawStyle( q , GLU_LINE );
    DrawSphere( p , 1.01*radi , n );
    gluQuadricDrawStyle( q , GLU_FILL );
  };

}

void Viewer::DrawBond(int atoma , int atomb , double radi , int n )
{
  double start[3],end[3],mid[3];

  start[0] = COORD[1][atoma][0];
  start[1] = COORD[1][atoma][1];
  start[2] = COORD[1][atoma][2];
  end[0] = COORD[1][atomb][0];
  end[1] = COORD[1][atomb][1];
  end[2] = COORD[1][atomb][2];
  mid[0] = (start[0]+end[0])/2.0;
  mid[1] = (start[1]+end[1])/2.0;
  mid[2] = (start[2]+end[2])/2.0;

  gluQuadricDrawStyle( q , GLU_FILL );
  glColor3dv( ELEMENT_COLOR[ATOMNO[atoma]] );
  DrawCylinder( start , mid , radi , n );

  glColor3dv( ELEMENT_COLOR[ATOMNO[atomb]] );
  DrawCylinder( mid, end , radi , n );
}


void Viewer::DrawMolecule( int ll, int ul, double srs, double crs)
{
  int iatom,jatom,ns;
  double dist,radi;
  double *posa,*posb,vec1[4],vec2[4];

  if (!drawmol) goto ARROW;

  ns = 50;
  glEnable( GL_LIGHTING );
  if ( srs == 0.0 && crs == 0.0 )
  {
    glDisable( GL_LIGHTING );
    glLineWidth( 2 );
    glBegin( GL_LINES );
  };

  for ( iatom = ll ; iatom <= ul ; iatom++ )
  {
    posa = &COORD[1][iatom][0]; 
    if ( srs != 0.0 )
    {
      if ( srs < 0.0 )
      {
        if ( srs == -0.0001 )
        {
          radi = ELEMENT_VDW_R[ATOMNO[iatom]];
        }
        else
        {
          radi = -srs;
        }
      }
      else
      {
        radi = srs*ELEMENT_COV_R[ATOMNO[iatom]];
      };
      if ( ATOMNO[iatom] == 0 ) 
      {
        radi = 0.10;
      };
      //cout << iatom << " " << radi << " " << ns << endl;
      DrawAtom( iatom , radi , ns );
    };

    if ( crs != 0.0 || srs == 0.0 )
    {
      //for ( jatom = iatom-1 ; jatom >= 0 ; jatom-- )
      for ( jatom = iatom-1 ; jatom >= ll ; jatom-- )
      {
        posb = &COORD[1][jatom][0]; 
        dist = ::distance(iatom,jatom);
        if ( dist < 1.3*(ELEMENT_COV_R[ATOMNO[iatom]] + ELEMENT_COV_R[ATOMNO[jatom]])&&
          ATOMNO[iatom] != 0 && ATOMNO[jatom] != 0 )
        {
          if ( srs == 0.0 && crs == 0.0 )
          {
            glColor3dv( ELEMENT_COLOR[ATOMNO[iatom]] );
            glVertex3d((double)posa[0],(double)posa[1],(double)posa[2]); 
            glColor3dv( ELEMENT_COLOR[ATOMNO[jatom]] );
            glVertex3d((double)posb[0],(double)posb[1],(double)posb[2]); 
          }
          else
          {
            if ( crs < 0.0 )
            {
              radi = - crs;
            }
            else
            {
              radi = ELEMENT_COV_R[ATOMNO[iatom]];
              if ( radi > ELEMENT_COV_R[ATOMNO[jatom]] )
              {
                radi = ELEMENT_COV_R[ATOMNO[jatom]];
              };
              radi *= crs;
            };
            DrawBond( iatom , jatom , radi , ns/2 ); 
          };
        };
      };
    };
  };
ARROW:
  DrawCPs();

  // Arrows
  glColor3f( arrow_color[0] , arrow_color[1] , arrow_color[2] );
  double s = 0.15;
  if (arrow_type == "Dipole") 
  {
    double norma = sqrt(pow(MOLDIPOLE[0],2) +
                        pow(MOLDIPOLE[1],2) +
                        pow(MOLDIPOLE[2],2));
    if (norma>1.0e-10) 
    {
        posa = &COORD[1][iatom][0]; 
      vec2[0] = MOLDIPOLE[0];
      vec2[1] = MOLDIPOLE[1];
      vec2[2] = MOLDIPOLE[2];
      vec1[0] = -vec2[0];
      vec1[1] = -vec2[1];
      vec1[2] = -vec2[2];
      DrawCylinder(vec1,vec2,s/2.0,10);
      vec1[0] = vec2[0] + s*MOLDIPOLE[0]/norma; 
      vec1[1] = vec2[1] + s*MOLDIPOLE[1]/norma; 
      vec1[2] = vec2[2] + s*MOLDIPOLE[2]/norma; 
      DrawCone(vec2,vec1,s,10);
    }
  }
  else if (arrow_type == "NFF HOMO" || arrow_type == "NFF LUMO" ||
           arrow_type == "NFFEX" || arrow_type == "NFFEY" ||
           arrow_type == "NFFEZ" || arrow_type == "Forces") 
  {
    double norma;
    double arrow[3];

    for ( iatom = ll ; iatom <= ul ; iatom++ )
    {
      posa = &COORD[1][iatom][0]; 
      if (arrow_type == "Forces") 
      {
        arrow[0] = 50.0*ATOMFORCE[iatom][0];
        arrow[1] = 50.0*ATOMFORCE[iatom][1];
        arrow[2] = 50.0*ATOMFORCE[iatom][2];
      }
      else if (arrow_type == "NFFEX") 
      {
        arrow[0] = ATOMNFFEX[iatom][0];
        arrow[1] = ATOMNFFEX[iatom][1];
        arrow[2] = ATOMNFFEX[iatom][2];
      }
      else if (arrow_type == "NFFEY") 
      {
        arrow[0] = ATOMNFFEY[iatom][0];
        arrow[1] = ATOMNFFEY[iatom][1];
        arrow[2] = ATOMNFFEY[iatom][2];
      }
      else if (arrow_type == "NFFEZ") 
      {
        arrow[0] = ATOMNFFEZ[iatom][0];
        arrow[1] = ATOMNFFEZ[iatom][1];
        arrow[2] = ATOMNFFEZ[iatom][2];
      }
      else if (arrow_type == "NFF HOMO") 
      {
        arrow[0] = ATOMNUCFUKUI_H[iatom][0];
        arrow[1] = ATOMNUCFUKUI_H[iatom][1];
        arrow[2] = ATOMNUCFUKUI_H[iatom][2];
      }
      else if (arrow_type == "NFF LUMO") 
      {
        arrow[0] = ATOMNUCFUKUI_L[iatom][0];
        arrow[1] = ATOMNUCFUKUI_L[iatom][1];
        arrow[2] = ATOMNUCFUKUI_L[iatom][2];
      }
      norma = sqrt(pow(arrow[0],2) + pow(arrow[1],2) + pow(arrow[2],2));
      if (norma>1.0e-10) 
      {
        posa = &COORD[1][iatom][0]; 
        vec1[0] = COORD[1][iatom][0] + 3.0*arrow[0]; 
        vec1[1] = COORD[1][iatom][1] + 3.0*arrow[1]; 
        vec1[2] = COORD[1][iatom][2] + 3.0*arrow[2]; 
        DrawCylinder(posa,vec1,s/2.0,10);
        vec2[0] = vec1[0] + s*arrow[0]/norma; 
        vec2[1] = vec1[1] + s*arrow[1]/norma; 
        vec2[2] = vec1[2] + s*arrow[2]/norma; 
        DrawCone(vec1,vec2,s,10);
      }
    }
  }

  if ( srs == 0.0 && crs == 0.0 )
  {
    glEnable( GL_LIGHTING );
    glLineWidth( 1 );
    glEnd();
  };
} 
void Viewer::DrawCPs( void )
{
  int i,ns;
  double radi;
  double p[3];

  ns = 20;
  radi = 0.12;

  glColor3f( cp_color[0],cp_color[1],cp_color[2]);
  for ( i = 0 ; i < NCPs ; i++ )
  {
    p[0] = CPCOORD[i][0];
    p[1] = CPCOORD[i][1];
    p[2] = CPCOORD[i][2];
    DrawSphere( p , radi , ns );
  };
} 

int Viewer::SelectedAtom( size_t n )
{
  for ( size_t i = 0 ; i < selected_atoms.size() ; i++ )
  {
    if ( n == selected_atoms[i] )
    {
      return i; 
    };
  };
  return -1; 
}

void Viewer::SaveSelection( int n )
{
  bool changed;

  changed = false;

  if ( prk->status & SFSURF )
  {
    prk->app->beep();
    prk->status &= ~SFSURF;
    QColor color; 
    color = QColorDialog::getColor( color , this );
    double col[4];
    pmgr->GetColor( n , col );
    col[0] = color.red()/255.0;
    col[1] = color.green()/255.0;
    col[2] = color.blue()/255.0;
    pmgr->SetColor( n , col );
    prk->SetCursor( "Normal" );
    Redraw();
  }
  else if ( prk->status & SFRV )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      prk->app->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 2 )
      {
        monitored_bonds.push_back( selected_atoms );
        selected_atoms.clear();
        prk->status &= ~SFRV;
        prk->Report("Atom %d selected, selection done",n + 1);  
        prk->SetCursor( "Normal" );
      }
      else
      {
        prk->Report("Atom %d selected, select another",n + 1);  
      };
      Redraw();
    };
  }

 else if ( prk->status & SFPG )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      prk->app->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 3 )
      {
        char gname[80];
        int ref[3];
        strcpy(gname,prk->geoeditor->gname);
        ref[0] = selected_atoms[0];
        ref[1] = selected_atoms[1];
        ref[2] = selected_atoms[2];
        prk->status &= ~SFPG;
        if (sgdrv(gname,ref))
        {
          prk->Report("Group building failed");
        }
        else
        {
          prk->geoeditor->Update();
          changed = true;
        }
      }
      else
      {
        Redraw();
      };
    };
  }
 else if ( prk->status & SFPB )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      prk->app->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 2 )
      {
        prk->status &= ~SFPB;
        changed = true;
      }
      else
      {
        prk->Report("Atom %d selected, select another",n + 1);
      };
      Redraw();
    };
  }
  else if ( prk->status & SFAV )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      prk->app->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 3 )
      {
        monitored_angles.push_back( selected_atoms );
        selected_atoms.clear();
        prk->status &= ~SFAV;
        prk->Report("Atom %d selected, selection done",n + 1);  
        prk->SetCursor( "Normal" );
      }
      else
      {
        prk->Report("Atom %d selected",n + 1);  
      };
      Redraw();
    };
  }
  else if ( prk->status & SFDV )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      prk->app->beep();
      selected_atoms.push_back( n );
      if ( selected_atoms.size() == 4 )
      {
        monitored_dihedrals.push_back( selected_atoms );
        selected_atoms.clear();
        prk->status &= ~SFDV;
        prk->Report("Atom %d selected, selection done",n + 1);  
        prk->SetCursor( "Normal" );
      }
      else
      {
        prk->Report("Atom %d selected",n + 1);  
      };
      Redraw();
    };
  } else   if ( prk->status & SFB )
  {
    prk->app->beep();
    prk->Report( "Atom %d selected",n+1);
    prk->status &= ~SFB;
    if ( NATOM == 2 )
    {
      prk->geoeditor->AddAtom( n + 1, 2 - n, 1);
      changed = true;
    }
    else
    {
      prk->status |= SFA;
      selected_atoms.push_back( n );
      Redraw();
      prk->Report("Atom %d selected",n + 1);
    };
  }
  else if ( prk->status & SFA )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      selected_atoms.push_back( n );
      prk->app->beep();
      prk->Report( "Atom %d selected",n+1);
      prk->status &= ~SFA;
      if ( NATOM == 3 )
      {
        prk->geoeditor->AddAtom( selected_atoms[0] + 1, n + 1,
                                   4 - selected_atoms[0] - n);
        changed = true;
      }
      else
      {
        prk->status |= SFD;
        Redraw();
        prk->Report( "Atom %d selected",n + 1);
      };
    };
  }
  else if ( prk->status & SFD )
  {
    if ( SelectedAtom( n ) < 0 )
    {
      selected_atoms.push_back( n );
      prk->app->beep();
      prk->Report("Atom %d selected",n+1);
      prk->status &= ~SFD;
      prk->geoeditor->AddAtom(selected_atoms[0] + 1,
                                selected_atoms[1] + 1, n + 1);
      changed = true;
    };
  };

  if ( changed )
  {
    selected_atoms.clear();
    prk->SetCursor( "Normal" );
  };
}

void Viewer::DrawSurfaces()
{
  for ( int i = 0 ; i < pmgr->NSurf(); i++ )
  {
    pmgr->PlotSurface( i );
  };
  // Extra call for drawing references
  pmgr->PlotSurface( (int) pmgr->NSurf() );
};

void Viewer::Write(char* text, Vector p, size_t fontsize)
{
  QFont font( "Helvetica" , fontsize );
  if ( bg[0] == 1.0 &&  bg[1] == 1.0 && bg[2] == 1.0 )
  {
    glColor3f( 0.0 , 0.0 , 0.0 );
  }
  else
  {
    glColor3f( 1.0 , 1.0 , 1.0 );
  }
  glDisable( GL_LIGHTING );
  glEnable( GL_LINE_STIPPLE );
  glLineStipple( 1 , 0xf0f0 );
  renderText( p[0] , p[1] , p[2] , QString(text), font );
  glDisable( GL_LINE_STIPPLE );
  glEnable( GL_LIGHTING );
}

void Viewer::DrawStrings()
{
  static int angle_step = 3;
  size_t i,j,iatom,jatom,katom;
  size_t fontsize=32;
  char word[MAX_STR_SIZE];
  double val,ss;
  Vector va,vb,vc,vd,vshift;
  Vector rva,rvb,rvc,rvd;
  Vector auxa,auxb,auxc;

  if ( bg[0] == 1.0 &&  bg[1] == 1.0 && bg[2] == 1.0 )
  {
   glColor3f( 0.0 , 0.0 , 0.0 );
  }
  else
  {
   glColor3f( 1.0 , 1.0 , 1.0 );
  }

  vshift = Vector( rot[2], rot[6] , rot[10] );

  glDisable( GL_LIGHTING );
  glEnable( GL_LINE_STIPPLE );
  glLineStipple( 1 , 0xf0f0 );

  strcpy( word , (image_render->label->text().toStdString()).c_str() );
  fontsize = 24;
  renderText( 10 , 50 , QString(word), QFont( "Helvetica" , fontsize ) );
  ss = 5.0*geocrs*scale;
  fontsize = size_t(MAX(4.0,18*scale));
  for ( i = 0 ; i < monitored_bonds.size() ; i++ )
  {
    iatom = monitored_bonds[i][0];
    jatom = monitored_bonds[i][1];
    va = Vector( &COORD[1][iatom][0] );
    vb = Vector( &COORD[1][jatom][0]);
    glBegin( GL_LINES );
      glVertex3d( va[0] , va[1] , va[2] );
      glVertex3d( vb[0] , vb[1] , vb[2] );
    glEnd();
    vc = va - vb;
    val = BohrToAngstrom(vc.Norm()); 
    va += vb;
    va *= 0.5;
    vb = vshift;
    vb *= ss;
    va += vb;
    sprintf(word,"%3.3f",val);
    Write(word, va, fontsize );
  };

  for ( i = 0 ; i < monitored_angles.size() ; i++ )
  {
    iatom = monitored_angles[i][0];
    jatom = monitored_angles[i][1];
    katom = monitored_angles[i][2];
    va = Vector(&COORD[1][iatom][0]);
    vb = Vector(&COORD[1][jatom][0]);
    vc = Vector(&COORD[1][katom][0]);

    rva = va - vb;
    rvc = vc - vb;
    rva ^= 1.0;
    rvc ^= 1.0;
    rvb = rva > rvc;
    val = 180.0/M_PI*acos( rva.Dot( rvc ) );
    sprintf(word,"%3.1f",val);
    rva ^= 1.0;
    rva.Rotate( rvb , val/2.0 );
    rvc = vb + rva;
    rva.Rotate( rvb , -val/2.0 );
    Write(word, rvc, fontsize );

    glBegin( GL_LINES );
    glVertex3d( va[0] , va[1], va[2] );
    glVertex3d( vb[0] , vb[1], vb[2] );
    glVertex3d( vb[0] , vb[1], vb[2] );
    glVertex3d( vc[0] , vc[1], vc[2] );
    glEnd();
    glBegin( GL_POINTS );
    for ( j = 0 ; j < (size_t) val ; j += angle_step )
    {
      rva.Rotate( rvb , (double)angle_step );
      rvc = vb + rva;
      glVertex3d( rvc[0] , rvc[1] , rvc[2] );
    };
    glEnd();
  };

  for ( i = 0 ; i < monitored_dihedrals.size() ; i++ )
  {
    iatom = monitored_dihedrals[i][0];
    va = Vector(&COORD[1][iatom][0]);
    iatom = monitored_dihedrals[i][1];
    vb = Vector(&COORD[1][iatom][0]);
    iatom = monitored_dihedrals[i][2];
    vc = Vector(&COORD[1][iatom][0]);
    iatom = monitored_dihedrals[i][3];
    vd = Vector(&COORD[1][iatom][0]);

    rva = va - vb;
    rvc = vc - vb;
    rvd = vd - vc;
    rva ^= 1.0;
    rvc ^= 1.0;
    rvd ^= 1.0;
    auxa = rvc;
    auxa *= -rva.Dot( rvc );
    auxa += rva;
    auxa ^= 1.0;
    auxb = rvc;
    auxb *= -rvd.Dot( rvc );
    auxb += rvd;
    auxb ^= 1.0;
    auxc = auxa;
    auxc.Rotate( rvc , -90.0 );
    val = auxa.Dot( auxb );
    if ( val < 0.0 )  
    {
      if ( auxb.Dot( auxc ) < 0.0 )
      {
        val = 360.0 - 180.0/M_PI*acos( val );
      }
      else
      {
        val =  180.0/M_PI*acos( val );
      }
    }
    else
    {
      if ( auxb.Dot( auxc ) < 0.0 )
      {
        val = 360.0 - 180.0/M_PI*acos( val );
      }
      else
      {
        val = 180.0/M_PI*acos( val );
      }
    };
    sprintf(word,"%3.1f",val);
    auxa = rvc;
    auxa *= -rva.Dot( rvc );
    auxa += rva;
    auxa.Rotate( rvc , -val/2.0 );
    auxb = vb + vc; 
    auxb *= 0.5; 
    auxb += auxa;
    auxa.Rotate( rvc , val/2.0 );
    Write(word, auxb,  fontsize);

    auxb = vb + vc; 
    auxb *= 0.5; 
    auxb += auxa;
    glBegin( GL_LINES );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( va[0] , va[1], va[2] );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vb[0] , vb[1], vb[2] );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vc[0] , vc[1], vc[2] );
    auxa.Rotate( rvc , -val );
    auxb = vb + vc; 
    auxb *= 0.5; 
    auxb += auxa;
    auxa.Rotate( rvc , val );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vb[0] , vb[1], vb[2] );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vc[0] , vc[1], vc[2] );
    glVertex3d( auxb[0] , auxb[1], auxb[2] );
    glVertex3d( vd[0] , vd[1], vd[2] );
    glEnd();

    glBegin( GL_POINTS );
    for ( j = 0 ; j < (size_t) val ; j += angle_step )
    {
      auxa.Rotate( rvc , -(double)angle_step );
      auxb = vb + vc; 
      auxb *= 0.5; 
      auxb += auxa;
      glVertex3d( auxb[0], auxb[1], auxb[2] );
    };
    glEnd();
  };

  for ( i = 0 ; i < user_label.size() ; i++ )
  {
    renderText( user_label[i].x, 
                user_label[i].y, 
                QString( user_label[i].text ), QFont( "Helvetica" , 32 )); 
  }

  // Atom labels are grey (Roberto Flores-Moreno, Feb 2008)
  glColor3f( tag_color[0] , tag_color[1] , tag_color[2] );
  if ( ! (tag_type == "None" ) )
  {
    ss = 10.0*geosrs*scale;
    vb = vshift;
    vb *= ss;
    for ( iatom = 0 ; iatom < (unsigned)NATOM ; iatom++ )
    {
      va = Vector(&COORD[1][iatom][0]);
      va += vb;
      if ( tag_type == "Symbols" )
      {
        sprintf(word,"%s",ELEMENT_SYMBOL[ATOMNO[iatom]]);
      }
      else if ( tag_type == "Numbers" )
      {
        sprintf(word,"(%d)",iatom+1);
      }
      else if ( tag_type == "Sym+Num" )
      {
        sprintf(word,"%s(%d)",ELEMENT_SYMBOL[ATOMNO[0]],iatom+1);
      }
      else if ( tag_type == "Charges" )
      {  
        sprintf(word,"%.2f",ATOMCHARGE[iatom]);
      }
      else if ( tag_type == "Fukui HOMO" )
      {  
        sprintf(word,"%.3f",ATOMFUKUI_H[iatom]);
      }
      else if ( tag_type == "Fukui LUMO" )
      {  
        sprintf(word,"%.3f",ATOMFUKUI_L[iatom]);
      }
      else if ( tag_type == "Fukui Average" )
      {  
        sprintf(word,"%.3f",(ATOMFUKUI_H[iatom]+ATOMFUKUI_L[iatom])/2.0);
      }
      else if ( tag_type == "Fukui Difference" )
      {  
        sprintf(word,"%.2f",(ATOMFUKUI_L[iatom]-ATOMFUKUI_H[iatom]));
      }
      else
      {
        break;
      }
      Write( word,va, fontsize );
    }
  };

  glEnable( GL_LIGHTING );
  glDisable( GL_LINE_STIPPLE );
}

void Viewer::ChangeBackgroundColor()
{
  QColor color; 
  color = QColorDialog::getColor( color , this );
  makeCurrent();
  bg[0] = color.red()/255.0;
  bg[1] = color.green()/255.0;
  bg[2] = color.blue()/255.0;
  glClearColor( bg[0] , bg[1] , bg[2] , bg[3] );
  Redraw();
}

void Viewer::ChangeTags( const QString &new_tag_type )
{
  tag_type = new_tag_type;
  Redraw();
}

void Viewer::ChangeArrows( const QString &new_arrow_type )
{
  arrow_type = new_arrow_type;
  Redraw();
}

void Viewer::ChangeTagColor()
{
  QColor color; 
  color = QColorDialog::getColor( color , this );
  tag_color[0] = color.red()/255.0;
  tag_color[1] = color.green()/255.0;
  tag_color[2] = color.blue()/255.0;
  Redraw();
}

void Viewer::ChangeArrowColor()
{
  QColor color; 
  color = QColorDialog::getColor( color , this );
  arrow_color[0] = color.red()/255.0;
  arrow_color[1] = color.green()/255.0;
  arrow_color[2] = color.blue()/255.0;
  Redraw();
}

void Viewer::DrawCone(double *start, double *end, 
                          double radi , int n )
{
  GLdouble v[3];
  GLdouble dval, ca, sa, cb, sb;
  GLdouble tmat[4][4];
  GLfloat height; 

  v[0] = (GLdouble)(end[0]-start[0]); 
  v[1] = (GLdouble)(end[1]-start[1]); 
  v[2] = (GLdouble)(end[2]-start[2]); 
  height = (GLfloat) sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);

  dval = sqrt(v[0]*v[0] + v[1]*v[1]); 
  if (dval != 0.0) 
  {
    ca = v[1]/dval;
    sa = v[0]/dval;
    v[1] = v[0]*sa + v[1]*ca;
  } 
  else
  {
    ca = 1.0;
    sa = 0.0;
  };

  dval = sqrt(v[1]*v[1] + v[2]*v[2]); 
  if (dval != 0.0) 
  {
    cb = v[2]/dval;
    sb = v[1]/dval;
  }
  else 
  {
    cb = 1.0;
    sb = 0.0;
  };

  tmat[0][0] = ca;   
  tmat[1][0] = sa*cb; 
  tmat[2][0] = sa*sb; 
  tmat[3][0] = (GLdouble)start[0];
  tmat[0][1] = (-sa); 
  tmat[1][1] = ca*cb; 
  tmat[2][1] = ca*sb; 
  tmat[3][1] = (GLdouble)start[1];
  tmat[0][2] = 0.0;   
  tmat[1][2] = (-sb); 
  tmat[2][2] = cb;    
  tmat[3][2] = (GLdouble)start[2];
  tmat[0][3] = 0.0;   
  tmat[1][3] = 0.0;   
  tmat[2][3] = 0.0;   
  tmat[3][3] = 1.0;

  glPushMatrix();
  glMultMatrixd((GLdouble*)tmat);
  gluCylinder( q, (GLfloat)radi, (GLfloat)0.0, height, n, 1);
  glPopMatrix();
}

